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

#include "StringBasics.h"
#include "Parameters.h"
#include <map>

struct StoredInfo
{
    int chrom;
    int start;
    int end;
    int totalReads;
    int numDups;
    int numQCFail;
    int numMapped;
    int numPaired;
    int numProper;
    int numZeroMapQ;
    int numLT10MapQ;
    int num255MapQ;
    int numMapQPass;
    double sumMapQ;
    int avgMapQCount;
    int depth;
    int numQ20;
};
    

bool readNextLine(IFILE inputFile, StoredInfo& nextLine, int& minChrom, int& minPos);
void updateSummary(const StoredInfo& nextLine, StoredInfo& sumLine);
bool writeSummary(IFILE outputFile, StoredInfo& summaryLine);
void initStoredInfo(StoredInfo& info);

void setupChromMap(std::map <std::string, int> &chromMap);

const char* chromNames[] = 
    {"1",
     "2",
     "3",
     "4",
     "5",
     "6",
     "7",
     "8",
     "9",
     "10",
     "11",
     "12",
     "13",
     "14",
     "15",
     "16",
     "17",
     "18",
     "19",
     "20",
     "21",
     "22",
     "X",
     "Y",
     "MT",
     "GL000207.1",
     "GL000226.1",
     "GL000229.1",
     "GL000231.1",
     "GL000210.1",
     "GL000239.1",
     "GL000235.1",
     "GL000201.1",
     "GL000247.1",
     "GL000245.1",
     "GL000197.1",
     "GL000203.1",
     "GL000246.1",
     "GL000249.1",
     "GL000196.1",
     "GL000248.1",
     "GL000244.1",
     "GL000238.1",
     "GL000202.1",
     "GL000234.1",
     "GL000232.1",
     "GL000206.1",
     "GL000240.1",
     "GL000236.1",
     "GL000241.1",
     "GL000243.1",
     "GL000242.1",
     "GL000230.1",
     "GL000237.1",
     "GL000233.1",
     "GL000204.1",
     "GL000198.1",
     "GL000208.1",
     "GL000191.1",
     "GL000227.1",
     "GL000228.1",
     "GL000214.1",
     "GL000221.1",
     "GL000209.1",
     "GL000218.1",
     "GL000220.1",
     "GL000213.1",
     "GL000211.1",
     "GL000199.1",
     "GL000217.1",
     "GL000216.1",
     "GL000215.1",
     "GL000205.1",
     "GL000219.1",
     "GL000224.1",
     "GL000223.1",
     "GL000195.1",
     "GL000212.1",
     "GL000222.1",
     "GL000200.1",
     "GL000193.1",
     "GL000194.1",
     "GL000225.1",
     "GL000192.1"
    };

std::map <std::string, int> chromMap;

std::map<std::string,int>::iterator chromMapIter;

char* chromBuffer = new char[100];

void usage()
{
    std::cerr << "Merge baseQC Cout-Based Summary Statistics.\n";
    std::cerr << "Usage: mergeStats --out <outputStatsFile> <inputStatsFiles>\n"
              << "\t--out output merged stats file\n"
              << "\tinputStatsFiles space separated list of files to merge.\n"
              << "\n";

}


int main(int argc, char ** argv)
{
    String output = "";
    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("out", &output)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));
    
    int numArgsProcessed = inputParameters.ReadWithTrailer(argc, argv);

    // numArgsProcessed does not include the program name, so add one to it.
    ++numArgsProcessed;

    if(numArgsProcessed == argc)
    {
        // No stats files specified
        std::cerr << "No stats files specified.\n";
        usage();
        return(-1);
    }

    if(output.Length() == 0)
    {
        std::cerr << "No output file sepcified, exiting.\n";
        usage();
        return(-1);
    }

    setupChromMap(chromMap);

    // Open the output file.
    IFILE outputFile = ifopen(output, "w");

    int numFiles = argc - numArgsProcessed;
    std::vector<IFILE> inputFiles;
    std::vector<StoredInfo> nextLine;
    inputFiles.resize(numFiles);
    nextLine.resize(numFiles);

    String header;
    String dataLine;

    int minChrom = 0x7FFFFFFF;
    int minPos = 0x7FFFFFFF;
    int nextMinChrom = 0x7FFFFFFF;
    int nextMinPos = 0x7FFFFFFF;

    for(int i = 0; i < numFiles; i++)
    {
        // Open the stats input files.
        inputFiles[i] = ifopen(argv[numArgsProcessed+i], "r");
        if(inputFiles[i] == NULL)
        {
            std::cerr << "Failed to open " << argv[numArgsProcessed+i] << " for reading.\n";
            exit(-1);
        }
        // Read the first line (this is the header).
        header.ReadLine(inputFiles[i]);

        // Read the first data line
        if(!readNextLine(inputFiles[i], nextLine[i], minChrom, minPos))
        {
            // there is not another record, so set chrom/pos to 0x7FFFFFFF
            nextLine[i].chrom = 0x7FFFFFFF;
            nextLine[i].start = 0x7FFFFFFF;
        }
    }

    // write the header.
    ifprintf(outputFile, "%s\n", header.c_str());

    StoredInfo sumLine;
    initStoredInfo(sumLine);

    bool done = false;
    while(!done)
    {
        done = true;
        // Keep looping through, merging records.
        for(int i = 0; i < numFiles; i++)
        {
            if(nextLine[i].start == 0x7FFFFFFF)
            {
                // This file is done processing, continue to the next one.
                continue;
            }
            done = false;

            // Check to see if this record matchs the min.
            if((minChrom == nextLine[i].chrom) && (minPos == nextLine[i].start))
            {
                // This is a min line, so accumulate
                updateSummary(nextLine[i], sumLine);
                // Used this line, so read the next line.
                if(!readNextLine(inputFiles[i], nextLine[i], nextMinChrom, nextMinPos))
                {
                    // there is not another record, so set chrom/pos to 0x7FFFFFFF
                    nextLine[i].chrom = 0x7FFFFFFF;
                    nextLine[i].start = 0x7FFFFFFF;
                }
            }
            // Check if it is the next min position.
            else if(nextLine[i].chrom < nextMinChrom)
            {
                // Lower Chromosome
                nextMinChrom = nextLine[i].chrom;
                nextMinPos = nextLine[i].start;
            }
            else if((nextLine[i].chrom == nextMinChrom) && (nextLine[i].start < nextMinPos))
            {
                // Lower position.
                nextMinPos = nextLine[i].start;
            }
        }
        if(!done)
        {
            writeSummary(outputFile, sumLine);
        }
        minChrom = nextMinChrom;
        minPos = nextMinPos;
        nextMinChrom = 0x7FFFFFFF;
        nextMinPos = 0x7FFFFFFF;
    }

    ifclose(outputFile);
    for(int i = 0; i < numFiles; i++)
    {
        ifclose(inputFiles[i]);
    }

    std::cerr << "Done writing to " << output << std::endl;

    return(0);
}


bool readNextLine(IFILE inputFile, StoredInfo& nextLine, int& minChrom, int& minPos)
{
    double avgMapQ = 0;

    static String dataLine;


    // Read the first data line
    if(dataLine.ReadLine(inputFile) < 0)
    {
        return(false);
    }
    // Parse the data line.
    if(sscanf(dataLine.c_str(), "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%d\t%d\t%d",
              chromBuffer, &(nextLine.start), &(nextLine.end), 
              &(nextLine.totalReads), &(nextLine.numDups), 
              &(nextLine.numQCFail), &(nextLine.numMapped),
              &(nextLine.numPaired), &(nextLine.numProper),
              &(nextLine.numZeroMapQ), &(nextLine.numLT10MapQ),
              &(nextLine.num255MapQ), &(nextLine.numMapQPass),
              &avgMapQ, &(nextLine.avgMapQCount), 
              &(nextLine.depth), &(nextLine.numQ20)) != 17)
    {
        std::cerr << "Failed reading line from " << inputFile->getFileName() << "\n";
        exit(-1);
    }
    
    // Convert the chromosome to it's integer value.
    chromMapIter = chromMap.find(chromBuffer);
    if(chromMapIter == chromMap.end())
    {
        std::cerr << "Skipping chromosome " << chromBuffer << std::endl;
        return(readNextLine(inputFile, nextLine, minChrom, minPos));
    }
    nextLine.chrom = chromMapIter->second;

    // Calculate the values for this data line.
    nextLine.sumMapQ = avgMapQ * nextLine.avgMapQCount;
    
    if(nextLine.chrom < minChrom)
    {
        // Lower Chromosome
        minChrom = nextLine.chrom;
        minPos = nextLine.start;
    }
    else if((nextLine.chrom == minChrom) && (nextLine.start < minPos))
    {
        // Lower position.
        minPos = nextLine.start;
    }
    return(true);
}

void updateSummary(const StoredInfo& nextLine, StoredInfo& sumLine)
{
    sumLine.chrom = nextLine.chrom;
    sumLine.start = nextLine.start;
    sumLine.end = nextLine.end;

    sumLine.totalReads += nextLine.totalReads;
    sumLine.numDups += nextLine.numDups;
    sumLine.numQCFail += nextLine.numQCFail;
    sumLine.numMapped += nextLine.numMapped;
    sumLine.numPaired += nextLine.numPaired;
    sumLine.numProper += nextLine.numProper;
    sumLine.numZeroMapQ += nextLine.numZeroMapQ;
    sumLine.numLT10MapQ += nextLine.numLT10MapQ;
    sumLine.num255MapQ += nextLine.num255MapQ;
    sumLine.numMapQPass += nextLine.numMapQPass;
    sumLine.sumMapQ += nextLine.sumMapQ;
    sumLine.avgMapQCount += nextLine.avgMapQCount;
    sumLine.depth += nextLine.depth;
    sumLine.numQ20 += nextLine.numQ20;
}


bool writeSummary(IFILE outputFile, StoredInfo& summaryLine)
{
    double avgMapQ = 0;

    if(summaryLine.avgMapQCount != 0)
    {
        avgMapQ = (double)(summaryLine.sumMapQ)/summaryLine.avgMapQCount;
    }
    ifprintf(outputFile,
             "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%d\t%d\t%d\n",
             chromNames[summaryLine.chrom], summaryLine.start, summaryLine.end,
             summaryLine.totalReads, summaryLine.numDups,
             summaryLine.numQCFail, summaryLine.numMapped,
             summaryLine.numPaired, summaryLine.numProper,
             summaryLine.numZeroMapQ, summaryLine.numLT10MapQ,
             summaryLine.num255MapQ, summaryLine.numMapQPass,
             avgMapQ, summaryLine.avgMapQCount,
             summaryLine.depth, summaryLine.numQ20);

    initStoredInfo(summaryLine);
    return(true);
}

void initStoredInfo(StoredInfo& info)
{
    info.chrom = 0;
    info.start = 0;
    info.end = 0;

    info.totalReads = 0;
    info.numDups = 0;
    info.numQCFail = 0;
    info.numMapped = 0;
    info.numPaired = 0;
    info.numProper = 0;
    info.numZeroMapQ = 0;
    info.numLT10MapQ = 0;
    info.num255MapQ = 0;
    info.numMapQPass = 0;
    info.sumMapQ = 0;
    info.avgMapQCount = 0;
    info.depth = 0;
    info.numQ20 = 0;
}


void setupChromMap(std::map <std::string, int> &chromMap)
{
    int mapIndex = 0;
    chromMap["1"] = mapIndex; ++mapIndex;
    chromMap["2"] = mapIndex; ++mapIndex;
    chromMap["3"] = mapIndex; ++mapIndex;
    chromMap["4"] = mapIndex; ++mapIndex;
    chromMap["5"] = mapIndex; ++mapIndex;
    chromMap["6"] = mapIndex; ++mapIndex;
    chromMap["7"] = mapIndex; ++mapIndex;
    chromMap["8"] = mapIndex; ++mapIndex;
    chromMap["9"] = mapIndex; ++mapIndex;
    chromMap["10"] = mapIndex; ++mapIndex;
    chromMap["11"] = mapIndex; ++mapIndex;
    chromMap["12"] = mapIndex; ++mapIndex;
    chromMap["13"] = mapIndex; ++mapIndex;
    chromMap["14"] = mapIndex; ++mapIndex;
    chromMap["15"] = mapIndex; ++mapIndex;
    chromMap["16"] = mapIndex; ++mapIndex;
    chromMap["17"] = mapIndex; ++mapIndex;
    chromMap["18"] = mapIndex; ++mapIndex;
    chromMap["19"] = mapIndex; ++mapIndex;
    chromMap["20"] = mapIndex; ++mapIndex;
    chromMap["21"] = mapIndex; ++mapIndex;
    chromMap["22"] = mapIndex; ++mapIndex;
    chromMap["X"] = mapIndex; ++mapIndex;
    chromMap["Y"] = mapIndex; ++mapIndex;
    chromMap["MT"] = mapIndex; ++mapIndex;
    chromMap["GL000207.1"] = mapIndex; ++mapIndex;
    chromMap["GL000226.1"] = mapIndex; ++mapIndex;
    chromMap["GL000229.1"] = mapIndex; ++mapIndex;
    chromMap["GL000231.1"] = mapIndex; ++mapIndex;
    chromMap["GL000210.1"] = mapIndex; ++mapIndex;
    chromMap["GL000239.1"] = mapIndex; ++mapIndex;
    chromMap["GL000235.1"] = mapIndex; ++mapIndex;
    chromMap["GL000201.1"] = mapIndex; ++mapIndex;
    chromMap["GL000247.1"] = mapIndex; ++mapIndex;
    chromMap["GL000245.1"] = mapIndex; ++mapIndex;
    chromMap["GL000197.1"] = mapIndex; ++mapIndex;
    chromMap["GL000203.1"] = mapIndex; ++mapIndex;
    chromMap["GL000246.1"] = mapIndex; ++mapIndex;
    chromMap["GL000249.1"] = mapIndex; ++mapIndex;
    chromMap["GL000196.1"] = mapIndex; ++mapIndex;
    chromMap["GL000248.1"] = mapIndex; ++mapIndex;
    chromMap["GL000244.1"] = mapIndex; ++mapIndex;
    chromMap["GL000238.1"] = mapIndex; ++mapIndex;
    chromMap["GL000202.1"] = mapIndex; ++mapIndex;
    chromMap["GL000234.1"] = mapIndex; ++mapIndex;
    chromMap["GL000232.1"] = mapIndex; ++mapIndex;
    chromMap["GL000206.1"] = mapIndex; ++mapIndex;
    chromMap["GL000240.1"] = mapIndex; ++mapIndex;
    chromMap["GL000236.1"] = mapIndex; ++mapIndex;
    chromMap["GL000241.1"] = mapIndex; ++mapIndex;
    chromMap["GL000243.1"] = mapIndex; ++mapIndex;
    chromMap["GL000242.1"] = mapIndex; ++mapIndex;
    chromMap["GL000230.1"] = mapIndex; ++mapIndex;
    chromMap["GL000237.1"] = mapIndex; ++mapIndex;
    chromMap["GL000233.1"] = mapIndex; ++mapIndex;
    chromMap["GL000204.1"] = mapIndex; ++mapIndex;
    chromMap["GL000198.1"] = mapIndex; ++mapIndex;
    chromMap["GL000208.1"] = mapIndex; ++mapIndex;
    chromMap["GL000191.1"] = mapIndex; ++mapIndex;
    chromMap["GL000227.1"] = mapIndex; ++mapIndex;
    chromMap["GL000228.1"] = mapIndex; ++mapIndex;
    chromMap["GL000214.1"] = mapIndex; ++mapIndex;
    chromMap["GL000221.1"] = mapIndex; ++mapIndex;
    chromMap["GL000209.1"] = mapIndex; ++mapIndex;
    chromMap["GL000218.1"] = mapIndex; ++mapIndex;
    chromMap["GL000220.1"] = mapIndex; ++mapIndex;
    chromMap["GL000213.1"] = mapIndex; ++mapIndex;
    chromMap["GL000211.1"] = mapIndex; ++mapIndex;
    chromMap["GL000199.1"] = mapIndex; ++mapIndex;
    chromMap["GL000217.1"] = mapIndex; ++mapIndex;
    chromMap["GL000216.1"] = mapIndex; ++mapIndex;
    chromMap["GL000215.1"] = mapIndex; ++mapIndex;
    chromMap["GL000205.1"] = mapIndex; ++mapIndex;
    chromMap["GL000219.1"] = mapIndex; ++mapIndex;
    chromMap["GL000224.1"] = mapIndex; ++mapIndex;
    chromMap["GL000223.1"] = mapIndex; ++mapIndex;
    chromMap["GL000195.1"] = mapIndex; ++mapIndex;
    chromMap["GL000212.1"] = mapIndex; ++mapIndex;
    chromMap["GL000222.1"] = mapIndex; ++mapIndex;
    chromMap["GL000200.1"] = mapIndex; ++mapIndex;
    chromMap["GL000193.1"] = mapIndex; ++mapIndex;
    chromMap["GL000194.1"] = mapIndex; ++mapIndex;
    chromMap["GL000225.1"] = mapIndex; ++mapIndex;
    chromMap["GL000192.1"] = mapIndex; ++mapIndex;
}

