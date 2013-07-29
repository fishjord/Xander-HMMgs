/*
 * Copyright (C) 2012 Jordan Fish <fishjord at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package edu.msu.cme.rdp.graph.cli;

import edu.msu.cme.rdp.alignment.hmm.HMMER3bParser;
import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import edu.msu.cme.rdp.graph.filter.BloomFilter;
import edu.msu.cme.rdp.graph.filter.BloomFilter.GraphState;
import edu.msu.cme.rdp.graph.search.HMMGraphSearch;
import edu.msu.cme.rdp.graph.search.SearchResult;
import edu.msu.cme.rdp.graph.search.SearchTarget;
import edu.msu.cme.rdp.graph.utils.BackTranslationIterator;
import edu.msu.cme.rdp.kmer.trie.ModelPositionKmerGenerator;
import edu.msu.cme.rdp.kmer.set.KmerIterator;
import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author fishjord
 */
public class BasicSearch {

    private final ProfileHMM forHMM, revHMM;
    private final BloomFilter bloom;
    private final HMMGraphSearch searcher;
    private final GraphState walker;

    public BasicSearch(final int maxPaths, final ProfileHMM forHMM, final ProfileHMM revHMM, final BloomFilter bloom) {
        this.forHMM = forHMM;
        this.revHMM = revHMM;
        this.searcher = new HMMGraphSearch(maxPaths);
        this.bloom = bloom;
        this.walker = bloom.new GraphState();
    }

    private List<SearchResult> processOneKmer(String queryId, char[] protKmer, int modelPos) throws InterruptedException {
        byte[] nuclMer;
        BackTranslationIterator it = new BackTranslationIterator(protKmer);
        List<SearchResult> ret = new ArrayList();
        int considered = 0;
        int explored = 0;
	long startTime = System.currentTimeMillis();

        while ((nuclMer = it.next()) != null) {
            walker.setStateDangerously(nuclMer);
            considered++;
            //if (!walker.hasCurrent()) {
            //    continue;
            //}

            explored++;

            /*SearchTarget target = new SearchTarget(forHMM.getName(),
                    queryId, "?", new String(nuclMer), 0,
                    modelPos - 1, forHMM, revHMM, bloom);

	    ret.addAll(searcher.search(target));*/
        }

	//System.err.println(queryId + "\t" + new String(protKmer) + "\t" + modelPos + "\t" + explored + "\t" + considered + "\t" + (System.currentTimeMillis() - startTime));

        return ret;
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 5) {
            System.err.println("USAGE: BasicSearch <bloom_filter> <for_hmm> <rev_hmm> <ref_file>");
            System.exit(1);
        }

        int k = Integer.valueOf(args[0]);

        File bloomFile = new File(args[1]);
        File forHMMFile = new File(args[2]);
        File revHMMFile = new File(args[3]);
        File refFile = new File(args[4]);

        File nuclOutFile = new File(refFile.getName() + "_nucl.fasta");
        File alignOutFile = new File(refFile.getName() + "_aln.fasta");
        File protOutFile = new File(refFile.getName() + "_prot.fasta");

        ProfileHMM forHMM;
        ProfileHMM revHMM;

        forHMM = HMMER3bParser.readModel(forHMMFile);
        revHMM = HMMER3bParser.readModel(revHMMFile);


        FastaWriter nuclOut = new FastaWriter(nuclOutFile);
        FastaWriter alignOut = new FastaWriter(alignOutFile);
        FastaWriter protOut = null;
        boolean isProt = forHMM.getAlphabet() == SequenceType.Protein;

        if (isProt) {
            protOut = new FastaWriter(protOutFile);
        }

        int contigCount = 1;

        long startTime;

        startTime = System.currentTimeMillis();
        ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(bloomFile)));
        BloomFilter bloom = (BloomFilter) ois.readObject();
        ois.close();

        BasicSearch search = new BasicSearch(k, forHMM, revHMM, bloom);

        System.err.println("Bloom filter loaded in " + (System.currentTimeMillis() - startTime) + " ms");

        System.err.println("Starting hmmgs search at " + new Date());
        System.err.println("*  Reference File:          " + refFile);
        System.err.println("*  Bloom file:              " + bloomFile);
        System.err.println("*  Forward hmm file:        " + forHMMFile);
        System.err.println("*  Reverse hmm file:        " + revHMMFile);
        System.err.println("*  Searching prot?:         " + isProt);
        System.err.println("*  # paths:                 " + k);
        System.err.println("*  Nucl contigs out file    " + nuclOutFile);
        System.err.println("*  Prot contigs out file    " + protOutFile);

        HMMBloomSearch.printHeader(System.out, isProt);

        Set<Kmer> processed = new HashSet();

        SequenceReader reader = new SequenceReader(refFile);
        Sequence seq;
        ModelPositionKmerGenerator iter;
        int ksize = bloom.getKmerSize();
	Set<String> seenProtMers = new HashSet();
	String s;

        try {
            while ((seq = reader.readNextSequence()) != null) {
		if(seq.getSeqName().startsWith("#")) {
		    continue;
		}
		long seq_time = System.currentTimeMillis();
                iter = new ModelPositionKmerGenerator(seq.getSeqString(), bloom.getKmerSize() / 3, (isProt ? SequenceType.Protein : SequenceType.Nucleotide));
                for (char[] protMer : iter) {
		    s = new String(protMer);
		    if(seenProtMers.contains(s)) {
			continue;
		    }
		    seenProtMers.add(s);
		    long t = System.currentTimeMillis();
                    for (SearchResult result : search.processOneKmer(seq.getSeqName(), protMer, iter.getModelPosition())) {
                        if(result.getNuclSeq().length() < ksize * 2) {
                            continue;
                        }
                        String seqid = "contig_" + (contigCount++);

                        HMMBloomSearch.printResult(seqid, isProt, result, System.out);

                        nuclOut.writeSeq(seqid, result.getNuclSeq());
                        alignOut.writeSeq(seqid, result.getAlignSeq());
                        if (isProt) {
                            protOut.writeSeq(seqid, result.getProtSeq());
                        }
                    }
                }
		System.err.println("Processed seq " + seq.getSeqName() + " in " + (System.currentTimeMillis() - seq_time) + " ms");
            }

        } finally {
            nuclOut.close();
            if (isProt) {
                protOut.close();
            }
            System.out.close();
        }

    }
}
