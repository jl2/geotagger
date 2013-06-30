/*
  main.cpp
  
  Copyright (c) 2013, Jeremiah LaRocco jeremiah.larocco@gmail.com

  Permission to use, copy, modify, and/or distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

/*
  Console application to geotag a directory of Jpeg files using a GPX file.
*/

#include <exiv2/exiv2.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <algorithm>

#include <ctime>

#define BOOST_SYSTEM_NO_DEPRECATED

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

#include "tinyxml2.h"
namespace tx = tinyxml2;


void findJpegFiles(fs::path directory, std::vector<fs::path> &jpegFiles) {
    // If the directory doesn't exist, it has no jpeg files, so return
    if (!fs::exists(directory)) return;

    std::copy_if(fs::directory_iterator(directory), fs::directory_iterator(),
                 std::back_inserter(jpegFiles), [](fs::path pn) {
                     std::string fname = pn.string();
                     bool rval = (fs::is_regular_file(pn) &&
                                  (boost::algorithm::ends_with(fname, ".jpg")
                                   || boost::algorithm::ends_with(fname, ".JPG")
                                   || boost::algorithm::ends_with(fname, ".jpeg")
                                   || boost::algorithm::ends_with(fname, ".JPEG")));
                     return rval;
                 });
}

void convertToDDMMSS(double val, int &degrees, int &minutes, double &seconds) {
    degrees = 0;
    minutes = 0;
    seconds = 0.0;

    double tempf = 0.0;
    if (val<0) {
        val = -val;
    }
    double rem =  std::modf(val, &tempf);
    degrees = int(tempf);
    if (rem < 0) {
        rem = -rem;
    }
    // std::cout << "degrees = " << degrees << " rem = " << rem << " rem*60 = " << (rem*60) << "\n";
    
    rem *= 60;
    double rem2 = std::modf(rem, &tempf);
    minutes = int(tempf);
    // std::cout << "minutes = " << minutes << " rem2 = " << rem2 << " rem2*60 = " << (rem2*60) << "\n";
    
    seconds = rem2 * 60;
}

struct GPXPoint{
    double lat;
    double lon;
    double ele;
    std::string dateStamp;
    int hour;
    int minute;
    int second;
    std::time_t timestamp;
};

void readGPX(std::string inFile, std::vector<GPXPoint> &out) {
    tx::XMLDocument doc;
    doc.LoadFile( inFile.c_str() );
    tx::XMLElement *gpx = doc.FirstChildElement("gpx");
    for (auto *trks = gpx->FirstChildElement("trk");
         trks;
         trks = trks->NextSiblingElement()) {
        for (auto segs = trks->FirstChildElement("trkseg");
             segs; segs = segs->NextSiblingElement()) {
            for (auto pt = segs->FirstChildElement("trkpt");
                 pt; pt = pt->NextSiblingElement()) {
                std::string lat = pt->Attribute("lat");
                std::string lon = pt->Attribute("lon");
                tx::XMLElement *elevationTag = pt->FirstChildElement("ele");
                tx::XMLElement *timeTag = pt->FirstChildElement("time");

                GPXPoint tmp;
                tmp.lat = std::stod(lat);
                tmp.lon = std::stod(lon);
                tmp.ele = std::stod(elevationTag->GetText());
                std::string fullTime = timeTag->GetText();
                tmp.dateStamp = fullTime.substr(0,4) + ":" + fullTime.substr(5,2) + ":" + fullTime.substr(8,2);
                std::tm timeData;
                // 2013-06-27T00:50:39Z
                strptime(fullTime.c_str(), "%Y-%m-%dT%H:%M:%SZ", &timeData);

                tmp.hour = timeData.tm_hour;
                tmp.minute = timeData.tm_min;
                tmp.second = timeData.tm_sec;
                // tmp.timestamp = std::mktime(&timeData);
                tmp.timestamp = timegm(&timeData);
                // std::cout <<  tmp.hour << ":" << tmp.minute << ":" << tmp.second << "\n";
                // std::cout << tmp.timestamp << "\n";
                out.push_back(tmp);
            }
        }
    }
    std::sort(out.begin(), out.end(),
              [](const GPXPoint &a,const GPXPoint &b) {
                  return a.timestamp < b.timestamp;
              });
}

bool findClosest(const std::vector<GPXPoint> &pts, std::time_t imgTime, size_t &idx) {
    idx = 0;

    size_t upper = pts.size();
    size_t lower = 0;
    std::tm *withTZone = std::gmtime(&imgTime);
    // std::cout << "Image: " << imgTime << " " << withTZone->tm_hour << ":" << withTZone->tm_min << ":" << withTZone->tm_sec << ")\n";

    withTZone = std::gmtime(&pts[0].timestamp);
    // std::cout << "Lower: " << pts[0].timestamp << " " << withTZone->tm_hour << ":" << withTZone->tm_min << ":" << withTZone->tm_sec << ")\n";
    withTZone = std::gmtime(&pts[upper-1].timestamp);
    // std::cout << "Upper: " << pts[upper-1].timestamp << " " << withTZone->tm_hour << ":" << withTZone->tm_min << ":" << withTZone->tm_sec << ")\n";

    if (imgTime < pts[0].timestamp) {
        // std::cout << imgTime << " < " <<  pts[0].timestamp<< "\n";
        return false;
        
    }
    if (imgTime > pts[upper-1].timestamp) {
        // std::cout << imgTime << " > " << pts[upper-1].timestamp << "\n";
        return false;
    }

    size_t curGuess = (upper-lower)/2 + lower;

    while (curGuess != upper && curGuess != lower) {
        // std::cout << "upper: " << upper << " guess: " << curGuess << " lower: " << lower << "\n";
        if (pts[curGuess].timestamp == imgTime) {
            idx = curGuess;
            return true;
        } else if (imgTime < pts[curGuess].timestamp ) {
            upper = curGuess;
            curGuess = (upper-lower)/2 + lower;
        } else if (imgTime > pts[curGuess].timestamp) {
            lower = curGuess;
            curGuess = (upper-lower)/2 + lower;
        }
    }

    idx = curGuess;
    return true;
}
std::time_t getImageTimeStamp(std::string tmpTime) {
    std::tm timeData;
    // 2013:06:26 16:33:44
    strptime(tmpTime.c_str(), "%Y:%m:%d %H:%M:%S", &timeData);
    auto tstamp = std::mktime(&timeData);
    std::tm *withTZone = std::localtime(&tstamp);
    tstamp = std::mktime(withTZone);
    return tstamp;
}

void clearGPSFields(Exiv2::ExifData &exifData) {
    std::list<std::string> fields = {"Exif.GPSInfo.GPSLatitude",
                                     "Exif.GPSInfo.GPSLongitude",
                                     "Exif.GPSInfo.GPSMapDatum",
                                     "Exif.GPSInfo.GPSLatitudeRef",
                                     "Exif.GPSInfo.GPSLongitudeRef",
                                     "Exif.GPSInfo.GPSLatitude",
                                     "Exif.GPSInfo.GPSLongitude",
                                     
                                     "Exif.GPSInfo.GPSAltitude",
                                     "Exif.GPSInfo.GPSAltitudeRef",

                                     "Exif.GPSInfo.GPSTimeStamp",
                                     "Exif.GPSInfo.GPSDateStamp",};
    for (auto fname : fields) {
        auto key = Exiv2::ExifKey(fname.c_str());
        auto iter = exifData.findKey(key);
        while (iter != exifData.end()) {
            exifData.erase(iter);
            iter = exifData.findKey(key);
        }
    }
}

int main(int argc, char* const argv[]) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <gpsfile> <directory>\n";
        return 1;
    }

    std::vector<GPXPoint> gpxData;
    readGPX(argv[1], gpxData);
    std::vector<fs::path> jpegFiles;

    findJpegFiles(argv[2], jpegFiles);

    for (const auto &fname : jpegFiles) {
        std::cout << "Processing: " << fname << "\n";
        try {
            Exiv2::Image::AutoPtr image = Exiv2::ImageFactory::open(fname.c_str());
            if (image.get() == 0) continue;
            image->readMetadata();

            Exiv2::ExifData &exifData = image->exifData();
            if (exifData.empty()) {
                std::string error = fname.string();
                error += ": No Exif data found in the file";
                throw Exiv2::Error(1, error);
            }
            std::string tmpTime = exifData["Exif.Photo.DateTimeOriginal"].toString();
            auto tstamp = getImageTimeStamp(tmpTime);

            size_t idx = 0;
            bool foundIt = findClosest(gpxData, tstamp, idx);

            if (!foundIt) {
                std::cout << fname << " is not on the GPX track!\n";
                continue;
            } else {
                std::cout << fname << " was at ("  << std::setprecision(10) << gpxData[idx].lat << ", " << std::setprecision(10) << gpxData[idx].lon << ")\n";
            }
            clearGPSFields(exifData);

            exifData["Exif.GPSInfo.GPSMapDatum"] = "WGS-84";

            exifData["Exif.GPSInfo.GPSAltitude"] = Exiv2::Rational(gpxData[idx].ele * 1000, 1000);
            exifData["Exif.GPSInfo.GPSAltitudeRef"] = Exiv2::byte(0);

            // Convert the latitude to DDD*MM'SS.SSS" and set
            int dd, mm;
            double ss;
            convertToDDMMSS(gpxData[idx].lat, dd, mm, ss);
            if (gpxData[idx].lat<0) {
                exifData["Exif.GPSInfo.GPSLatitudeRef"] = "S";
            } else {
                exifData["Exif.GPSInfo.GPSLatitudeRef"] = "N";
            }

            Exiv2::URationalValue::AutoPtr latitude(new Exiv2::URationalValue);
            latitude->value_.push_back(std::make_pair(dd,1));
            latitude->value_.push_back(std::make_pair(mm,1));
            latitude->value_.push_back(std::make_pair(std::trunc(ss*10000)-1,10000));
            auto latKey = Exiv2::ExifKey("Exif.GPSInfo.GPSLatitude");
            exifData.add(latKey, latitude.get());

            convertToDDMMSS(gpxData[idx].lon, dd, mm, ss);
            Exiv2::URationalValue::AutoPtr longitude(new Exiv2::URationalValue);
            if (gpxData[idx].lon<0) {
                exifData["Exif.GPSInfo.GPSLongitudeRef"] = "W";
            } else {
                exifData["Exif.GPSInfo.GPSLongitudeRef"] = "E";
            }
            longitude->value_.push_back(std::make_pair(dd,1));
            longitude->value_.push_back(std::make_pair(mm,1));
            longitude->value_.push_back(std::make_pair(int(ss*10000)-1,10000));
            auto longKey = Exiv2::ExifKey("Exif.GPSInfo.GPSLongitude");
            exifData.add(longKey, longitude.get());


            Exiv2::URationalValue::AutoPtr timestamp(new Exiv2::URationalValue);
            timestamp->value_.push_back(std::make_pair(gpxData[idx].hour,1));
            timestamp->value_.push_back(std::make_pair(gpxData[idx].minute,1));
            timestamp->value_.push_back(std::make_pair(gpxData[idx].second,1));

            auto timeKey = Exiv2::ExifKey("Exif.GPSInfo.GPSTimeStamp");
            exifData.add(timeKey, timestamp.get());
            
            exifData["Exif.GPSInfo.GPSDateStamp"] = gpxData[idx].dateStamp.c_str();

            image->setExifData(exifData);
            image->writeMetadata();
        }
        catch (Exiv2::AnyError& e) {
            std::cout << "Caught Exiv2 exception '" << e.what() << "'\n";
            continue;
        }
    }
    return 0;
}
