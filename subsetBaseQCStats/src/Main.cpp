/*
 *  Copyright (C) 2011  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Parameters.h"
#include "NonOverlapRegions.h"

int readRegions(String& regions, NonOverlapRegions& regionList);
bool writeLine(IFILE outputFile);

const unsigned int BUFFER_SIZE = 1000;
char readBuffer[BUFFER_SIZE];

// Just for the chromosome buffer.
const unsigned int CHROM_BUFFER_SIZE = 100;

int main(int argc, char ** argv)
{
    String input;
    String regions;
    String output;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("inStats", &input)
        LONG_STRINGPARAMETER("regionList", &regions)
        LONG_STRINGPARAMETER("outStats", &output)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));
    
    inputParameters.Read(argc, argv);

    // Check for required parameters.
    if(input.IsEmpty() || regions.IsEmpty() || output.IsEmpty())
    {
        // The required parameters were not specified.
        std::cerr << "Narrow down the stats to just a subset of positions.\n";
        std::cerr << "Usage: subsetBaseQCStats --in <originalStatsFile> --regionList <subset of regions> --outStats <outputStatsFile>\n"
                  << "\n";
        std::cerr << "\t\t--inStats    : stats file to narrow down to just a subset of positions" << std::endl;
        std::cerr << "\t\t--regionList : File containing the subset of regions to keep (assumed to be sorted)\n"
                  << "\t\t               Formated as chr<tab>start_pos<tab>end_pos.\n" 
                  << "\t\t               Positions are 0 based and the end_pos is not included in the region." << std::endl;
        return(-1);
    }

    IFILE inStats = ifopen(input, "r");
    if(inStats == NULL)
    {
        std::cerr << "Failed to open input stats file: " << input << std::endl;
        return(-1);
    }
    IFILE outStats = ifopen(output, "w");
    if(outStats == NULL)
    {
        std::cerr << "Failed to open output stats file: " << output 
                  << std::endl;
        ifclose(inStats);
        return(-1);
    }

    NonOverlapRegions regionList;
    int regionStat = readRegions(regions, regionList);

    if(regionStat != 0)
    {
        ifclose(inStats);
        ifclose(outStats);
        return(regionStat);
    }
    
    // Files were successfully opened, so reading the input file.
    // Buffer for reading the input lines into.
    bool error = false;
    char chrom[CHROM_BUFFER_SIZE];
    int pos;
    bool firstLine = true;
    // Keep reading the input file until the end is reached.
    while(!inStats->ifgetline(readBuffer, BUFFER_SIZE))
    {
        // Read a line from the file, parsing it to get the position.
        if(sscanf(readBuffer, "%s\t%d", chrom, &pos) != 2)
        {
            // Failed to read the line.
            if(firstLine)
            {
                // Header line.
                error &= writeLine(outStats);
                firstLine = false;
            }
            else
            {
                std::cerr << "Failed to read one of the stats lines from the input file.\n";
                error = true;
            }
            continue;
        }
        // Successfully read/parsed the line, so check if it is in the region.
        if(regionList.inRegion(chrom, pos))
        {
            error &= writeLine(outStats);
        }
    }

    // Done reading the input file.
    ifclose(inStats);
    ifclose(outStats);

    if(error)
    {
        return(-1);
    }
    return(0);
}

 
int readRegions(String& regions, NonOverlapRegions& regionList)
{
    IFILE inRegions = ifopen(regions, "r");
    if(inRegions == NULL)
    {
        std::cerr << "Failed to open input regions file: " << regions
                  << std::endl;
        return(-1);
    }

    // File was successfully opened, so read the regions.
    while(!inRegions->ifgetline(readBuffer, BUFFER_SIZE))
    {
        char* chrom = strtok(readBuffer,"\t");
        char* startStr = strtok(NULL,"\t");
        char* endStr = strtok(NULL,"\t");
        if((chrom != NULL) && (startStr != NULL) && (endStr != NULL))
        {
            // Successfully read a line.
            regionList.add(chrom, atoi(startStr), atoi(endStr));
        }
        else
        {
            // Line not properly formatted.
            std::cerr << "Invalid Line found in region list, continuing.\n";
        }
    }
    return(0);
}


bool writeLine(IFILE outputFile)
{
    static unsigned int writeLen = 0;
    
    if(outputFile == NULL)
    {
        return(false);
    }

    // In the region, so write the buffer.
    writeLen = strlen(readBuffer);
    if(writeLen+1 < BUFFER_SIZE)
    {
        readBuffer[writeLen] = '\n';
        ++writeLen;
    }
    if(outputFile->ifwrite(readBuffer, writeLen) != writeLen)
    {
        std::cerr << "subsetBaseQCStats: Failed to write a line to the output file.\n";
        return(false);
    }
    return(true);
}
