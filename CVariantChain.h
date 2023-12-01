/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CVariantChain.h
 * Author: mwittig
 *
 * Created on July 22, 2019, 9:56 AM
 */

#ifndef CVARIANTCHAIN_H
#define CVARIANTCHAIN_H




class CVariantChain {
public:
    CVariantChain(int maxThreads=2): m_isbt_anno(NULL), m_maxThreads(maxThreads), m_activeThreads(0) {}
    CVariantChain(CISBTAnno* isbtAnno, int maxThreads=5): m_isbt_anno(isbtAnno), m_maxThreads(maxThreads), m_activeThreads(0){}
    CVariantChain(const CVariantChain& orig);
    CVariantChain& operator=(const CVariantChain& orig);
    virtual ~CVariantChain();
    
    friend std::ostream& operator<<(std::ostream& os, const CVariantChain& me);
    
    /// this function is usually used to add ISBT variant for genomic positions
    /// where LRG sequence is != reference sequence (e.g del261G in ABO as the del is the hg19 reference)
    /// these SNPs are always homozyguous for the alternative allele
    void addHA(const CIsbtVariant& var);
    //void addHR(const CIsbtVariant& var);
    bool add(const CVcfSnp& var);
    
    std::set<CIsbtGt> getPossibleGenotypes()const;
    
    // this is used in the test data generating module. do not call this function
    void removeHA(const CIsbtVariant& var);
    
    std::map<std::string,set<CVariantChainVariation>>& getChains(){return m_chains;}
    
    
private:
    
    // std::set<CIsbtGt>&, CIsbtGtAllele, CIsbtGtAllele,map<string,set<CVariantChainVariation>>::const_iterator, int
    void getPossibleGenotypesMT(std::set<CIsbtGt>& vars, CIsbtGtAllele allele_A, CIsbtGtAllele allele_B, 
            map<string,set<CVariantChainVariation>>::const_iterator iter, int type = 0)const;
   
    void runInThread(
        std::function<void(std::set<CIsbtGt>&, CIsbtGtAllele, CIsbtGtAllele,
                           std::map<std::string, std::set<CVariantChainVariation>>::const_iterator, int)> func,
        std::set<CIsbtGt>& vars, CIsbtGtAllele allele_A, CIsbtGtAllele allele_B,
        std::map<std::string, std::set<CVariantChainVariation>>::const_iterator iter,
        int type = 0) const { // Hinzugefügtes const für runInThread
        std::unique_lock<std::mutex> lock(m_mutex);

        // Warten, bis die Anzahl der aktiven Threads kleiner als m_maxThreads ist
        m_condition.wait(lock, [this] { return m_activeThreads < m_maxThreads; });

        // Erhöhen Sie die Anzahl der aktiven Threads
        ++m_activeThreads;

        // Starten Sie einen neuen Thread, um die Funktion auszuführen
        std::thread([this, func, &vars, allele_A, allele_B, iter, type]() {
            func(vars, allele_A, allele_B, iter, type);

            // Reduzieren Sie die Anzahl der aktiven Threads und benachrichtigen Sie andere Threads
            std::unique_lock<std::mutex> lock(m_mutex);
            --m_activeThreads;
            m_condition.notify_one();
        }).detach();
    }
    void waitForCompletion()const {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_condition.wait(lock, [this] { return m_activeThreads == 0; });
    }
    
    CISBTAnno* m_isbt_anno;
    mutable int m_maxThreads;
    mutable int m_activeThreads;
    mutable std::mutex m_mutex;
    mutable std::condition_variable m_condition;
    /// key is the phasing chain id with two presets
    /// ha == homozygous for the alternative
    /// hr == homozygous for the reference
    /// no == non phased
    /// everything else is a unique phasing identifier
    std::map<std::string,set<CVariantChainVariation>> m_chains;

    
};

#endif /* CVARIANTCHAIN_H */

