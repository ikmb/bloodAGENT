/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.cc to edit this template
 */

/* 
 * File:   CScoreHaplotype.cpp
 * Author: mwittig
 * 
 * Created on December 1, 2023, 3:28 PM
 */

#include <iostream>
#include <vector>
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <libgen.h>
#include <limits>

#include <regex>

#include "CIsbtVariant.h"
#include "CBigWigReader.h"
#include "CScoreHaplotype.h"

CScoreHaplotype::~CScoreHaplotype() {
}

int CScoreHaplotype::performAlignment(){
    const int typicalGapPenalty = -1;  // Typische Strafe für Lücken
    const int specialGapPenalty = -2;   // Höhere Strafe, wenn "A" involviert
    const int matchScore = 2;   // Punkte für Übereinstimmungen
    const int typicalMismatchPenalty = -1; // Typische Strafe für Nicht-Übereinstimmungen
    const int specialMismatchPenalty = -2; // Höhere Strafe, wenn "A" involviert

    vector<CIsbtVariant> haplotype1,haplotype2;
    haplotype1.assign(m_haplotype1.begin(), m_haplotype1.end());
    haplotype2.assign(m_haplotype2.begin(), m_haplotype2.end());
    const int rows = haplotype1.size() + 1;
    const int cols = haplotype2.size() + 1;

    // Initialisiere die Matrix für die dynamische Programmierung
    std::vector<std::vector<int>> dp(rows, std::vector<int>(cols, 0));

    int maxScore = 0;
    [[maybe_unused]] int maxI = 0, maxJ = 0;

    // Fülle die Matrix und finde den maximalen Score
    // Fülle die Matrix und finde den maximalen Score
    for (int i = 1; i < rows; ++i) {
        for (int j = 1; j < cols; ++j) {
            int gapPenalty = (haplotype1[i - 1].isHighImpactSNP() || haplotype2[j - 1].isHighImpactSNP()) ? specialGapPenalty : typicalGapPenalty;
            int mismatchPenalty = (haplotype1[i - 1].isHighImpactSNP() || haplotype2[j - 1].isHighImpactSNP()) ? specialMismatchPenalty : typicalMismatchPenalty;

            int scoreDiagonal = dp[i - 1][j - 1] + (haplotype1[i - 1] == haplotype2[j - 1] ? matchScore : mismatchPenalty);
            int scoreUp = dp[i - 1][j] + gapPenalty;
            int scoreLeft = dp[i][j - 1] + gapPenalty;

            dp[i][j] = std::max(0, std::max(scoreDiagonal, std::max(scoreUp, scoreLeft)));

            if (dp[i][j] > maxScore) {
                maxScore = dp[i][j];
                maxI = i;
                maxJ = j;
            }
        }
    }

    // Hier kannst du den maximalen Score oder die ausgerichteten Sequenzen weiterverarbeiten
    std::cerr << "Maximaler Score: " << maxScore << std::endl;

    /*
    // Beispiel: Ausgabe der ausgerichteten Sequenzen
    std::string alignedSequence1, alignedSequence2;

    while (maxI > 0 && maxJ > 0 && dp[maxI][maxJ] != 0) {
        if (dp[maxI][maxJ] == dp[maxI - 1][maxJ - 1] + (m_haplotype1[maxI - 1] == m_haplotype2[maxJ - 1] ? matchScore : mismatchPenalty)) {
            alignedSequence1 = m_haplotype1[maxI - 1] + alignedSequence1;
            alignedSequence2 = m_haplotype2[maxJ - 1] + alignedSequence2;
            --maxI;
            --maxJ;
        } else if (dp[maxI][maxJ] == dp[maxI - 1][maxJ] + gapPenalty) {
            alignedSequence1 = m_haplotype1[maxI - 1] + alignedSequence1;
            alignedSequence2 = '-' + alignedSequence2;
            --maxI;
        } else {
            alignedSequence1 = '-' + alignedSequence1;
            alignedSequence2 = m_haplotype2[maxJ - 1] + alignedSequence2;
            --maxJ;
        }
    }

    std::cout << "Aligned Sequence 1: " << alignedSequence1 << std::endl;
    std::cout << "Aligned Sequence 2: " << alignedSequence2 << std::endl;
    std::cout << "Alignment durchgeführt." << std::endl;
   */
    return maxScore;
}
