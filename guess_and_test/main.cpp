/*
Copyright(c) 2012, Ilya Vorobyev und Vasiliy Usatyuk
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met :
*Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and / or other materials provided with the distribution.
* Neither the name of the <organization> nor the
names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include".\myLib\regularLDPC.h"

vector<vector<int> > PROTOGRAPH;
//= {
//    { 0, 0, 0, 1, 0, 0, 0, 1, 2 },
//    { 0, 0, 0, 0, 0, 1, 1, 1, 2 },
//    { 0, 1, 0, 0, 1, 1, 1, 0, 2 },
//    { 0, 0, 1, 0, 1, 1, 0, 1, 2 },
//    { 2, 1, 1, 1, 1, 0, 1, 0, 1 }};
ll CIRCULANT_SIZE;
int VARIABLE_NODES;
int CHECK_NODES;
ll DESIRED_NUMBER_OF_MATRICES = -1;

vector<vector<vector<int> > > gen(int EIRA) {
    vector<vector<vector<int> > > a(CHECK_NODES, vector<vector<int> >(VARIABLE_NODES));
    if (EIRA == 0) {
        for (int r = 0; r < CHECK_NODES; ++r) {
            for (int c = 0; c < VARIABLE_NODES; ++c) {
                a[r][c].resize(PROTOGRAPH[r][c]);
                while (true) {
                    for (int k = 0; k < PROTOGRAPH[r][c]; ++k) {
                        a[r][c][k] = rand() % CIRCULANT_SIZE;
                    }
                    sort(a[r][c].begin(), a[r][c].end());
                    bool ok = 1;
                    for (size_t i = 1; i < a[r][c].size(); ++i) {
                        if (a[r][c][i] == a[r][c][i - 1]) {
                            ok = 0;
                            break;
                        }
                    }
                    if (ok)
                        break;
                }
            }
        }
        return a;
    }
    if (EIRA > 0) {//only 0/1 in protograph 
        for (int r = 0; r < CHECK_NODES; ++r) {
            for (int c = 0; c < VARIABLE_NODES; ++c) {
                a[r][c].resize(PROTOGRAPH[r][c]);
                if (PROTOGRAPH[r][c] == 0)
                    continue;
                if (CHECK_NODES + c >= VARIABLE_NODES + 1)
                    a[r][c][0] = 0;
                else
                    a[r][c][0] = rand() % CIRCULANT_SIZE;
            }
        }

        if (EIRA == 1) {
            a[CHECK_NODES - 1][VARIABLE_NODES - CHECK_NODES][0] = a[0][VARIABLE_NODES - CHECK_NODES][0];
            for (int i = 1; i + 1 < CHECK_NODES; ++i) {
                if (PROTOGRAPH[i][VARIABLE_NODES - CHECK_NODES] == 1)
                    a[i][VARIABLE_NODES - CHECK_NODES][0] = 0;
            }
        }
        if (EIRA == 2) {
            a[CHECK_NODES - 1][VARIABLE_NODES - CHECK_NODES][0] = a[0][VARIABLE_NODES - CHECK_NODES][0] = 0;
        }
        return a;
    }
}



struct entry {
    int r, c, id;
    int index, rowIndex, colIndex;
    entry() {}
    entry(int _r, int _c, int _id) : r(_r), c(_c), id(_id) {}
    entry(int _r, int _c, int _id, int _index, int _rowIndex, int _colIndex) : r(_r), c(_c), id(_id), index(_index), rowIndex(_rowIndex), colIndex(_colIndex) {}

};

vector<entry> entries;
vector<vector<entry> > entriesByRow, entriesByCol;
vector<vector<vector<entry> > > entriesByRowAndCol;

void enumerateEntries() {
    entriesByRow.resize(CHECK_NODES);
    entriesByCol.resize(VARIABLE_NODES);
    entriesByRowAndCol.resize(CHECK_NODES);
    for (int i = 0; i < CHECK_NODES; ++i)
        entriesByRowAndCol[i].resize(VARIABLE_NODES);
    for (int r = 0; r < CHECK_NODES; ++r) {
        for (int c = 0; c < VARIABLE_NODES; ++c) {
            for (int i = 0; i < PROTOGRAPH[r][c]; ++i) {
                entry x(r, c, i, entries.size(), entriesByRow[r].size(), entriesByCol[c].size());
                entries.push_back(x);
                entriesByRow[r].push_back(x);
                entriesByCol[c].push_back(x);
                entriesByRowAndCol[r][c].push_back(x);
            }
        }
    }
}


struct CycleEnum {
    vector<entry> cycle;
    int SZ;
    CycleEnum(int _SZ) {
        SZ = _SZ;
        cycle.resize(SZ);
    }
    bool next(int step = -1) {
        if (step == -1)
            step += SZ;
        if (step == SZ - 1) {
            int r = cycle[step - 1].r;
            int c = cycle[0].c;
            for (int i = cycle[step].id + 1; i < PROTOGRAPH[r][c]; ++i) {
                if ((c == cycle[step - 1].c) && (i == cycle[step - 1].id))
                    continue;
                if ((r == cycle[0].r) && (i == cycle[0].id))
                    continue;
                cycle[step] = entriesByRowAndCol[r][c][i];
                return 1;
            }
            return next(step - 1);
        }
        if (step == 0) {
            for (size_t i = cycle[0].index + 1; i < entries.size(); ++i) {
                cycle[0] = entries[i];
                if (init(1))
                    return 1;
            }
            return 0;
        }
        if (step & 1) {
            int r = cycle[step - 1].r;
            for (size_t i = cycle[step].rowIndex + 1; i < entriesByRow[r].size(); ++i) {
                if (cycle[step - 1].rowIndex == i)
                    continue;
                cycle[step] = entriesByRow[r][i];
                if (init(step + 1))
                    return 1;
            }
            return next(step - 1);
        }
        int c = cycle[step - 1].c;
        for (size_t i = cycle[step].colIndex + 1; i < entriesByCol[c].size(); ++i) {
            if (cycle[step - 1].colIndex == i)
                continue;
            cycle[step] = entriesByCol[c][i];
            if (init(step + 1))
                return 1;
        }
        return next(step - 1);
    }
    bool init(int step) {
        if (step == SZ)
            return 1;
        if (step + 1 == SZ) {
            int r = cycle[step - 1].r;
            int c = cycle[0].c;
            for (int i = 0; i < PROTOGRAPH[r][c]; ++i) {
                if ((c == cycle[step - 1].c) && (i == cycle[step - 1].id))
                    continue;
                if ((r == cycle[0].r) && (i == cycle[0].id))
                    continue;
                cycle[step] = entriesByRowAndCol[r][c][i];
                return 1;
            }
            return 0;
        }
        if (step & 1) {
            int r = cycle[step - 1].r;
            for (size_t i = 0; i < entriesByRow[r].size(); ++i) {
                if (i == cycle[step - 1].rowIndex)
                    continue;
                cycle[step] = entriesByRow[r][i];
                if (init(step + 1))
                    return 1;
            }
            return 0;
        }
        int c = cycle[step - 1].c;
        for (size_t i = 0; i < entriesByCol[c].size(); ++i) {
            if (i == cycle[step - 1].colIndex)
                continue;
            cycle[step] = entriesByCol[c][i];
            if (init(step + 1))
                return 1;
        }
        return 0;
    }
    bool init() {
        size_t id = 0;
        while (id < entries.size()) {
            cycle[0] = entries[id];
            if (init(1))
                break;
            ++id;
        }
        return id < entries.size();
    }
};

void print(const vector<vector<vector<int> > >& a) {
    for (int i = 0; i < CHECK_NODES; ++i) {
        for (int j = 0; j < VARIABLE_NODES; ++j) {
            for (size_t k = 0; k + 1 < a[i][j].size(); ++k)
                cout << a[i][j][k] << "&";
            if (a[i][j].empty())
                cout << -1 << "\t";
            else
                cout << a[i][j][a[i][j].size() - 1] << "\t";
        }
        cout << endl;
    }
}

void eprint(const vector<vector<vector<int> > >& a) {
    for (int i = 0; i < CHECK_NODES; ++i) {
        for (int j = 0; j < VARIABLE_NODES; ++j) {
            for (size_t k = 0; k + 1 < a[i][j].size(); ++k)
                cerr << a[i][j][k] << "&";
            if (a[i][j].empty())
                cerr << -1 << "\t";
            else
                cerr << a[i][j][a[i][j].size() - 1] << "\t";
        }
        cerr << endl;
    }
}

//bool girthAtLeast6Manual(const vector<vector<vector<int> > >& a) {
//    //cycles inside one entry
//    for (int i = 0; i < CHECK_NODES; ++i) {
//        for (int j = 0; j < VARIABLE_NODES; ++j) {
//            for (size_t i1 = 0; i1 < a[i][j].size(); ++i1) {
//                for (size_t i2 = i1 + 1; i2 < a[i][j].size(); ++i2) {
//                    for (size_t i3 = i1; i3 < a[i][j].size(); ++i3) {
//                        if (i3 == i2)
//                            continue;
//                        for (size_t i4 = i1 + 1; i4 < a[i][j].size(); ++i4) {
//                            if (i3 == i4)
//                                continue;
//                            if ((a[i][j][i1] + a[i][j][i3] - a[i][j][i2] - a[i][j][i4] + 2 * CIRCULANT_SIZE) % CIRCULANT_SIZE == 0)
//                                return 0;
//                        }
//                    }
//                }
//            }
//        }
//    }
//    //cycles between two entries in one row
//    for (int i = 0; i < CHECK_NODES; ++i) {
//        for (int i1 = 0; i1 < VARIABLE_NODES; ++i1) {
//            if (a[i][i1].size() < 2)
//                continue;
//            set<int> dif;
//            for (size_t i11 = 0; i11 < a[i][i1].size(); ++i11) {
//                for (size_t i12 = i11 + 1; i12 < a[i][i1].size(); ++i12) {
//                    dif.insert(a[i][i1][i12] - a[i][i1][i11]);
//                    dif.insert(CIRCULANT_SIZE + a[i][i1][i11] - a[i][i1][i12]);
//                }
//            }
//            for (int i2 = i1 + 1; i2 < VARIABLE_NODES; ++i2) {
//                if (a[i][i2].size() < 2)
//                    continue;
//                for (size_t i21 = 0; i21 < a[i][i2].size(); ++i21) {
//                    for (size_t i22 = i21 + 1; i22 < a[i][i2].size(); ++i22) {
//                        if (dif.find(a[i][i2][i22] - a[i][i2][i21]) != dif.end())
//                            return 0;
//                    }
//                }
//            }
//        }
//    }
//    //cycles between two entries in one column
//    for (int i = 0; i < VARIABLE_NODES; ++i) {
//        for (int i1 = 0; i1 < CHECK_NODES; ++i1) {
//            if (a[i1][i].size() < 2)
//                continue;
//            set<int> dif;
//            for (size_t i11 = 0; i11 < a[i1][i].size(); ++i11) {
//                for (size_t i12 = i11 + 1; i12 < a[i1][i].size(); ++i12) {
//                    dif.insert(a[i1][i][i12] - a[i1][i][i11]);
//                    dif.insert(CIRCULANT_SIZE + a[i1][i][i11] - a[i1][i][i12]);
//                }
//            }
//            for (int i2 = i1 + 1; i2 < CHECK_NODES; ++i2) {
//                if (a[i2][i].size() < 2)
//                    continue;
//                for (size_t i21 = 0; i21 < a[i2][i].size(); ++i21) {
//                    for (size_t i22 = i21 + 1; i22 < a[i2][i].size(); ++i22) {
//                        if (dif.find(a[i2][i][i22] - a[i2][i][i21]) != dif.end())
//                            return 0;
//                    }
//                }
//            }
//        }
//    }
//    //cycles between four entries
//    for (int r1 = 0; r1 < CHECK_NODES; ++r1) {
//        for (int c1 = 0; c1 < VARIABLE_NODES; ++c1) {
//            if (a[r1][c1].empty())
//                continue;
//            for (int r2 = r1 + 1; r2 < CHECK_NODES; ++r2) {
//                if (a[r2][c1].empty())
//                    continue;
//                for (int c2 = c1 + 1; c2 < VARIABLE_NODES; ++c2) {
//                    if ((a[r1][c2].empty()) || (a[r2][c2].empty()))
//                        continue;
//                    for (size_t i11 = 0; i11 < a[r1][c1].size(); ++i11) {
//                        for (size_t i12 = 0; i12 < a[r1][c2].size(); ++i12) {
//                            for (size_t i21 = 0; i21 < a[r2][c1].size(); ++i21) {
//                                for (size_t i22 = 0; i22 < a[r2][c2].size(); ++i22) {
//                                    if ((a[r1][c1][i11] + a[r2][c2][i22] - a[r1][c2][i12] - a[r2][c1][i21] + 2 * CIRCULANT_SIZE) % CIRCULANT_SIZE == 0)
//                                        return 0;
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    return 1;
//}

bool girthAtLeast6Enum(const vector<vector<vector<int> > >& a) {
    CycleEnum cycle(4);
    if (!cycle.init())
        return 1;
    while (true) {
        int res = 0;
        for (int i = 0; i * 2 < 4; ++i) {
            res = res + CIRCULANT_SIZE + a[cycle.cycle[2 * i].r][cycle.cycle[2 * i].c][cycle.cycle[2 * i].id] - a[cycle.cycle[2 * i + 1].r][cycle.cycle[2 * i + 1].c][cycle.cycle[2 * i + 1].id];
        }
        if (res % CIRCULANT_SIZE == 0) {
            /*for (int i = 0; i < 4; ++i) {
            cout << cycle.cycle[i].r << " " << cycle.cycle[i].c << " " << cycle.cycle[i].id << endl;
            }*/
            return 0;

        }
        if (!cycle.next())
            return 1;
    }
}

//void test() {
//    time_t start = time(NULL);
//    ll iterationCount = 0;
//    while (true) {
//        ++iterationCount;
//        vector<vector<vector<int> > > a = gen();
//        bool res1 = girthAtLeast6Manual(a);
//        bool res2 = girthAtLeast6Enum(a);
//        if (res1 != res2) {
//            cerr << "ERROR" << endl;
//            cout << res1 << endl;
//            cout << res2 << endl;
//            print(a);
//            break;
//        }
//        if (iterationCount % 1000000 == 0) {
//            cout << iterationCount / 1000000 << " million iterations\n";
//            cout << time(NULL) - start << " seconds\n";
//        }
//    }
//
//}

bool noCycles(int g, const vector<vector<vector<int> > >& a) {
    CycleEnum cycle(g);
    if (!cycle.init())
        return 1;
    while (true) {
        ll res = 0;
        for (int i = 0; i * 2 < g; ++i) {
            res = res + CIRCULANT_SIZE + a[cycle.cycle[2 * i].r][cycle.cycle[2 * i].c][cycle.cycle[2 * i].id] - a[cycle.cycle[2 * i + 1].r][cycle.cycle[2 * i + 1].c][cycle.cycle[2 * i + 1].id];
        }
        if (res % CIRCULANT_SIZE == 0) {
            /*for (int i = 0; i < 4; ++i) {
            cout << cycle.cycle[i].r << " " << cycle.cycle[i].c << " " << cycle.cycle[i].id << endl;
            }*/
            return 0;

        }
        if (!cycle.next())
            return 1;
    }

}

int getGirth(const vector<vector<vector<int> > >& a) {
    for (int g = 4; g < 1000; g += 2) {
        if (!noCycles(g, a))
            return g;
    }
    cerr << "girth >= 1000\n";
    return (1 << 31) - 1;
}





int main(int argc, char* argv[]) {
    bool validInput = 1;
    //if (argc != 11) {
    //    validInput = 0;
    //    //std::cerr << "Usage: " << argv[0] << " -seed SEED -girth GIRTH -circulant CIRCULANT_SIZE -file FILENAME" << std::endl;
    //    //return 1;

    //}
    ll SEED = -1;
    ll GIRTH = -1;
    //ll CIRCULANT_SIZE = -1;
    int EIRA = 0;//0-no, 1-EIRA,2-EIRA_B

    string INPUT_FILENAME = "";
    for (int i = 1; i + 1 < argc; ++i) {
        if (string(argv[i]) == "-seed") {
            validInput = validInput && toUnsignedInt(argv[i + 1], SEED);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-EIRA") {
            EIRA = 1;
            continue;
        }
        if (string(argv[i]) == "-EIRA_B") {
            EIRA = 2;
            continue;
        }
        if (string(argv[i]) == "-girth") {
            validInput = validInput && toUnsignedInt(argv[i + 1], GIRTH);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-circulant") {
            validInput = validInput && toUnsignedInt(argv[i + 1], CIRCULANT_SIZE);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-numberOfMatrices") {
            validInput = validInput && toUnsignedInt(argv[i + 1], DESIRED_NUMBER_OF_MATRICES);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-file") {
            INPUT_FILENAME = argv[i + 1];
            ++i;
            continue;
        }

    }
    if ((GIRTH < 0) || (SEED < 0) || (CIRCULANT_SIZE < 0) || (DESIRED_NUMBER_OF_MATRICES < 0) || (INPUT_FILENAME == ""))
        validInput = 0;
    if (!validInput) {
        std::cerr << "Usage: " << argv[0] << " -seed SEED -girth GIRTH -circulant CIRCULANT_SIZE -numberOfMatrices DESIRED_NUMBER_OF_MATRICES -file INPUT_FILENAME [-EIRA or -EIRA_B]" << std::endl;
        return 1;
    }
    if (GIRTH & 1) {
        cerr << "girth must be even\n";
        return 1;
    }
    srand(SEED);
    freopen(INPUT_FILENAME.c_str(), "r", stdin);

    cin >> VARIABLE_NODES >> CHECK_NODES;
    PROTOGRAPH.assign(CHECK_NODES, vector<int>(VARIABLE_NODES));
    for (int i = 0; i < CHECK_NODES; ++i) {
        for (int j = 0; j < VARIABLE_NODES; ++j) {
            cin >> PROTOGRAPH[i][j];
        }
    }
    fclose(stdin);
    string folderName = toStr(CHECK_NODES) + "_" + toStr(VARIABLE_NODES) + "_" + toStr(CIRCULANT_SIZE) + "girth" + toStr(GIRTH);
    system(("mkdir " + folderName).c_str());
    string outputFilename = folderName + "/" + toStr(CHECK_NODES) + "_" + toStr(VARIABLE_NODES) + "_" + toStr(CIRCULANT_SIZE) + "girth" + toStr(GIRTH) + "seed" + toStr(SEED) + 
        "protograph_from_" + INPUT_FILENAME + "_matrix";
    enumerateEntries();

    time_t start = time(NULL);
    ll iterationCount = 0;
    //int res = 0;
    //vector<vector<vector<int> > > mtr;
    ll successCount = 0;
    while (successCount < DESIRED_NUMBER_OF_MATRICES) {
        ++iterationCount;
        vector<vector<vector<int> > > a = gen(EIRA);
        int g = getGirth(a);
        if (g >= GIRTH) {
            freopen((outputFilename + toStr(iterationCount) + ".txt").c_str(), "w", stdout);
            ++successCount;
            //res = g;
            //mtr = a;
            //cout << "girth = " << g << endl;
            cout << VARIABLE_NODES << "\t" << CHECK_NODES << "\t" << CIRCULANT_SIZE << endl;
            print(a);
            //cout << iterationCount << " iterations\n";
            //cout << time(NULL) - start << " seconds\n";
            cout << endl;
            fclose(stdout);
            cerr << "girth = " << g << endl;
            eprint(a);
            cerr << iterationCount << " iterations\n";
            cerr << time(NULL) - start << " seconds\n";
            cerr << endl;           
        }
        if (iterationCount % 1000000000 == 0) {
            cerr << iterationCount / 1000000000 << " billions iterations\n";
            cerr << time(NULL) - start << " seconds\n";
            cerr << endl;
        }
    }
    return 0;
}