/*
 * File:   PgSATests.h
 * Author: Tomek
 *
 * Created on 6 luty 2014, 18:19
 */

#ifndef PGSATESTS_H
#define	PGSATESTS_H

#include "testdata.h"
#include "../index/PgSAIndexFactory.h"

using namespace PgSAReadsSet;
using namespace PgSAIndex;

namespace PgSATest {

    int repeat = 11;
    int testKmersCount = 100000;

    int kValue;
    
    PgSAIndexStandard* pgsaIndex;
    string* testkmers = 0;
    std::pair<t_reads_c, t_read>* testFactors;
    

    typedef DefaultPgSAIndex<uint_reads_cnt_std, unsigned int
        , uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std,
            DefaultSuffixArrayOfConstantLengthTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, 4>::Type> PgSAIndexStandardImpl;

//    void drawDataK(ConcatenatedReadsSet* readsSet, int k) {
//        testkmers = generateTextFactors(readsSet, k, testKmersCount);
//    }

    void doFilterTTTs() {
        int tttsCount = 0;
        for (int i = 0; i < testKmersCount; i++) {
            if (testkmers[i].find_first_of("ACGN") == string::npos) {
                tttsCount++;
                int idx = (i + testKmersCount - 1) % testKmersCount;
                testkmers[i] = testkmers[idx];
                testFactors[i] = testFactors[idx];
            }
        }
        cout << "Filtered " << tttsCount << " TTTT...TTT.\n";
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void drawDataK(ReadsSetInterface<uint_read_len, uint_reads_cnt>* readsSet, int k) {
        if (testkmers != 0)
            delete[]testkmers;
        if (testFactors != 0)
            delete[]testFactors;
            
        testFactors = generateTagFactors(readsSet->readsCountVirtual(), readsSet->maxReadLengthVirtual(), k, testKmersCount);
        testkmers = generateTextFactors(readsSet, k, testFactors, testKmersCount);
        
        kValue = k;
//        cout << testkmers[testKmersCount - 1] << "\n";

    }

    inline const string getKmer(t_reads_c i) {
        const std::pair<t_reads_c, t_read> pair = testFactors[i];
        return pgsaIndex->getReadVirtual(pair.first).substr(pair.second, kValue);
    }
    
    template<typename uint_read_len, typename uint_reads_cnt>
    void drawDataK(ReadsSetInterface<uint_read_len, uint_reads_cnt>* readsSet, int k, bool pairScrambled ) {
        if (testkmers != 0)
            delete[]testkmers;
        if (testFactors != 0)
            delete[]testFactors;
            
        if (pairScrambled) 
            testFactors = generatePairScrambledTagFactors(readsSet->readsCountVirtual(), readsSet->maxReadLengthVirtual(), k, testKmersCount);
        else        
            testFactors = generateTagFactors(readsSet->readsCountVirtual(), readsSet->maxReadLengthVirtual(), k, testKmersCount); 
            
        testkmers = generateTextFactors(readsSet, k, testFactors, testKmersCount);
        
        kValue = k;   
    }

    void printResult(string testid, uint_max occurrencesCount, vector<int> times, uint_max res) {
        std::sort(times.begin(), times.end());
        int time_min = times[0];
        int time_med = times[times.size() / 2];

        if (occurrencesCount)
            printf("%-4s%-3d%-13s%10s%14s%11s%14s\n", testid.substr(0,3).c_str(), (int) (testkmers[0].length()), toString(res).c_str(),
                   toString(time_med * 1000000.0L / testKmersCount, 0).c_str(),
                   toString(time_med * 1000000.0L / occurrencesCount, 3).c_str(),
                   toString(time_min * 1000000.0L / testKmersCount, 0).c_str(),
                   toString(time_min * 1000000.0L / occurrencesCount, 3).c_str());
    }

    void q1Test(int repeat, vector<int>& times, uint_max &res) {
        times.clear();
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();
            vector<t_reads_c> q1res;

            t_reads_c occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                pgsaIndex->reportReads(testkmers[i], q1res);
                occurrencesCount += q1res.size();
            }

            times.push_back(clock_millis());
            res = occurrencesCount;
        }
    }

    void q2Test(int repeat, vector<int>& times, uint_max &res) {
        times.clear();
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();

            t_reads_c occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                occurrencesCount += pgsaIndex->countReads(testkmers[i]);
            }
            times.push_back(clock_millis());
            res = occurrencesCount;
        }

    }

    /*
                int tmp = inHowManyReadsDoesFOccur(testkmers[i]);
                occurrencesCount += tmp;
                cout << "Found " << testkmers[i] << " in " << tmp << " reads\n";
    */

    void q3Test(int repeat, vector<int>& times, uint_max &res) {
        times.clear();
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();
            vector<StandardOccurrence> q3res;

            uint_max occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                pgsaIndex->reportOccurrences(testkmers[i], q3res);
                occurrencesCount += q3res.size();
            }
            times.push_back(clock_millis());
            res = occurrencesCount;
        }

    }

    void q4Test(int repeat, vector<int>& times, uint_max &res) {
        times.clear();
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();

            uint_max occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                occurrencesCount += pgsaIndex->countOccurrences(testkmers[i]);
                cout << testkmers[i] << ": " << occurrencesCount << "\n";
            }

            times.push_back(clock_millis());
            res = occurrencesCount;
        }
    }

    void q4TestImpl(int repeat, vector<int>& times, uint_max &res) {
        PgSAIndexStandardImpl* index = (PgSAIndexStandardImpl*) pgsaIndex;

        times.clear();
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();

            uint_max occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                occurrencesCount += index->countOccurrencesImpl(testkmers[i]);
//                cout << testkmers[i] << ": " << occurrencesCount << "\n";
            }

            times.push_back(clock_millis());
            res = occurrencesCount;
        }

    }

     void q5Test(int repeat, vector<int>& times, uint_max &res) {
        times.clear();
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();
            vector<t_reads_c> q1res;

            t_reads_c occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                pgsaIndex->reportReadsWithSingleOccurrence(testkmers[i], q1res);
                occurrencesCount += q1res.size();
            }

            times.push_back(clock_millis());
            res = occurrencesCount;
        }
    }

    void q6Test(int repeat, vector<int>& times, uint_max &res) {
        times.clear();
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();

            t_reads_c occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                occurrencesCount += pgsaIndex->countSingleOccurrences(testkmers[i]);
            }
            times.push_back(clock_millis());
            res = occurrencesCount;
        }

    }

    void q7Test(int repeat, vector<int>& times, uint_max &res) {
        times.clear();
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();
            vector<StandardOccurrence> q3res;

            uint_max occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                pgsaIndex->reportSingleOccurrences(testkmers[i], q3res);
                occurrencesCount += q3res.size();
            }
            times.push_back(clock_millis());
            res = occurrencesCount;
        }

    }
    
    void q1pTest(int repeat, vector<int>& times, uint_max &res) {
        times.clear();
        
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();
            vector<t_reads_c> q1res;

            t_reads_c occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                pgsaIndex->reportReads(testFactors[i].first, testFactors[i].second, kValue, q1res);
                occurrencesCount += q1res.size();
            }

            times.push_back(clock_millis());
            res = occurrencesCount;
        }
    }

    void q2pTest(int repeat, vector<int>& times, uint_max &res) {
        times.clear();
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();

            t_reads_c occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                occurrencesCount += pgsaIndex->countReads(testFactors[i].first, testFactors[i].second, kValue);
            }
            times.push_back(clock_millis());
            res = occurrencesCount;
        }

    }

    /*
                int tmp = inHowManyReadsDoesFOccur(testkmers[i]);
                occurrencesCount += tmp;
                cout << "Found " << testkmers[i] << " in " << tmp << " reads\n";
    */

    void q3pTest(int repeat, vector<int>& times, uint_max &res) {
        times.clear();
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();
            vector<StandardOccurrence> q3res;

            uint_max occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                pgsaIndex->reportOccurrences(testFactors[i].first, testFactors[i].second, kValue, q3res);
                occurrencesCount += q3res.size();
            }
            times.push_back(clock_millis());
            res = occurrencesCount;
        }

    }

    void q4pTest(int repeat, vector<int>& times, uint_max &res) {
        times.clear();
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();

            uint_max occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                occurrencesCount += pgsaIndex->countOccurrences(testFactors[i].first, testFactors[i].second, kValue);
//                cout << testkmers[i] << ": " << occurrencesCount << "\n";
            }

            times.push_back(clock_millis());
            res = occurrencesCount;
        }
    }

    void q4pTestImpl(int repeat, vector<int>& times, uint_max &res) {
        PgSAIndexStandardImpl* index = (PgSAIndexStandardImpl*) pgsaIndex;

        times.clear();
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();

            uint_max occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                occurrencesCount += index->countOccurrencesImpl(testFactors[i].first, testFactors[i].second, kValue);
//                cout << testkmers[i] << ": " << occurrencesCount << "\n";
            }

            times.push_back(clock_millis());
            res = occurrencesCount;
        }

    }

     void q5pTest(int repeat, vector<int>& times, uint_max &res) {
        times.clear();
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();
            vector<t_reads_c> q1res;

            t_reads_c occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                pgsaIndex->reportReadsWithSingleOccurrence(testFactors[i].first, testFactors[i].second, kValue, q1res);
                occurrencesCount += q1res.size();
            }

            times.push_back(clock_millis());
            res = occurrencesCount;
        }
    }

    void q6pTest(int repeat, vector<int>& times, uint_max &res) {
        times.clear();
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();

            t_reads_c occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                occurrencesCount += pgsaIndex->countSingleOccurrences(testFactors[i].first, testFactors[i].second, kValue);
            }
            times.push_back(clock_millis());
            res = occurrencesCount;
        }

    }

    void q7pTest(int repeat, vector<int>& times, uint_max &res) {
        times.clear();
        for(int j = 0; j < repeat; j++) {
            clock_checkpoint();
            vector<StandardOccurrence> q3res;

            uint_max occurrencesCount = 0;
            for(int i = 0; i < testKmersCount; i++) {
                pgsaIndex->reportSingleOccurrences(testFactors[i].first, testFactors[i].second, kValue, q3res);
                occurrencesCount += q3res.size();
            }
            times.push_back(clock_millis());
            res = occurrencesCount;
        }

    }
    
    void runTestQueries(int repeat, bool byPosition) {
        vector<int> times;
        uint_max occurrencesCount;
        uint_max res;
        q4Test(1, times, occurrencesCount);

        cout << "\n *** TEST STARTED *** \nTotal occurences |Q3|: " << occurrencesCount << "\n";
        cout << "Queries about " << testKmersCount << " randomly chosen k-mers repeated " << repeat << " times.\n\n";

        cout << "Q  k  resCount       time med  time/occ med   time min  time/occ min   [nsec]\n\n";

        if (!byPosition) {
/*        
            q1Test(repeat, times, res);
            printResult("Q1 (In which reads does f occur?)", occurrencesCount, times, res);

            q2Test(repeat, times, res);
            printResult("Q2 (In how many reads does f occur?)", occurrencesCount, times, res);

            q3Test(repeat, times, res);
            printResult("Q3 (What are the occurrence positions of f?)", occurrencesCount, times, res);
*/
            q4Test(repeat, times, res);
            printResult("Q4 (What is the number of occurrences of f?)", occurrencesCount, times, res);
/*
    //        q4TestImpl(repeat, times, res);
    //        printResult("Q4I (No virtual calls) (What is the number of occurrences of f?)", occurrencesCount, times, res);

            q5Test(repeat, times, res);
            printResult("Q5 (In which reads does f occur only once?)", occurrencesCount, times, res);

            q6Test(repeat, times, res);
            printResult("Q6 (In how many reads does f occur only once?)", occurrencesCount, times, res);

            q7Test(repeat, times, res);
            printResult("Q7 (What are the occurrence positions of f in the reads where it occurs only once?)", occurrencesCount, times, res);
*/

        } else {
        
            q1pTest(repeat, times, res);
            printResult("Q1p(In which reads does f occur?)", occurrencesCount, times, res);

            q2pTest(repeat, times, res);
            printResult("Q2p(In how many reads does f occur?)", occurrencesCount, times, res);

            q3pTest(repeat, times, res);
            printResult("Q3p(What are the occurrence positions of f?)", occurrencesCount, times, res);

            q4pTest(repeat, times, res);
            printResult("Q4p(What is the number of occurrences of f?)", occurrencesCount, times, res);

    //        q4pTestImpl(repeat, times, res);
    //        printResult("Q4Ip(No virtual calls) (What is the number of occurrences of f?)", occurrencesCount, times, res);

            q5pTest(repeat, times, res);
            printResult("Q5p(In which reads does f occur only once?)", occurrencesCount, times, res);

            q6pTest(repeat, times, res);
            printResult("Q6p(In how many reads does f occur only once?)", occurrencesCount, times, res);

            q7pTest(repeat, times, res);
            printResult("Q7p(What are the occurrence positions of f in the reads where it occurs only once?)", occurrencesCount, times, res);

        }
            
        cout << "\n\n";
    }

    void runTests(PgSAIndexStandard* index) {
        pgsaIndex = index;
        cout << "*****************************************************************************\n";
        cout << index->getDescription();
//        cout << "\n";

        drawDataK(index, 11);
        runTestQueries(repeat, false);

        drawDataK(index, 16);
        runTestQueries(repeat, false);

        drawDataK(index, 22);
        runTestQueries(repeat, false);

        delete[]testkmers;
        testkmers = 0;
    }

    void runTest(PgSAIndexStandard* index, int _repeat, int _testKmersCount, int k, bool pairScrambled, bool filterTTTs, bool byPosition) {
        pgsaIndex = index;
        repeat = _repeat;
        testKmersCount = _testKmersCount;
        drawDataK(index, k, pairScrambled);
        if (filterTTTs)
            doFilterTTTs();
        runTestQueries(repeat, byPosition);
        delete[]testkmers;
        testkmers = 0;
    }

    void runTests(PgSAIndexStandard* index, int _repeat, int _testKmersCount) {
        repeat = _repeat;
        testKmersCount = _testKmersCount;
        runTests(index);
    }
}

#endif	/* PGSATESTS_H */

