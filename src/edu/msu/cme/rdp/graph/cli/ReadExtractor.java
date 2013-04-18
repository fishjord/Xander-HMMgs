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

import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.kmer.set.KmerIterator;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import java.io.*;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 *
 * @author fishjord
 */
public class ReadExtractor {

    private final int k;
    private final Set<Kmer> kmerSet;
    private int processedSeqs = 0;
    private int writtenSeqs = 0;


    public ReadExtractor(File contigFile, int k) throws IOException {
        this.k = k;

        List<Sequence> contigSeqs = SequenceReader.readFully(contigFile);

        kmerSet = new HashSet();

        Sequence seq;
        for (int index = 0; index < contigSeqs.size(); index++) {
            seq = contigSeqs.get(index);
            String seqstr = seq.getSeqString();
            addKmers(seqstr, index);
        }
    }

    private void addKmers(String seqString, int contigIndex) {
        KmerIterator kmerGen = new KmerIterator(seqString, k);
        Kmer val;

        while (kmerGen.hasNext()) {
            val = kmerGen.next();
            kmerSet.add(val);
        }
    }

    private void processRead(Sequence seq, FastaWriter out) throws IOException {
        if(seq.getSeqString().length() < k) {
            return;
        }

        if(checkSeqForKnownKmers(seq.getSeqString())) {
            out.writeSeq(seq);
            writtenSeqs++;
        } else {
            String rc = IUBUtilities.reverseComplement(seq.getSeqString());
            if(checkSeqForKnownKmers(rc)) {
                out.writeSeq(seq.getSeqName() + "_rc", rc);
                writtenSeqs++;
            }
        }
        processedSeqs++;
    }

    private boolean checkSeqForKnownKmers(String seqString) {
        KmerIterator kmerGen = new KmerIterator(seqString, k);
        Kmer val;
        Set<String> kmers;

        while (kmerGen.hasNext()) {
            val = kmerGen.next();

            if (kmerSet.contains(val)) {
                return true;
            }
        }

        return false;
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 3 && args.length != 4) {
            System.err.println("USAGE: ReadKmerMapper <nucl_contig_file> <reads_file> <k> [#threads]");
            System.exit(1);
        }

        final File nuclContigs = new File(args[0]);
        final File readsFile = new File(args[1]);
        final int k = Integer.valueOf(args[2]);
        final int maxThreads;
        final int maxTasks = 25000;

        if (args.length == 4) {
            maxThreads = Integer.valueOf(args[3]);
        } else {
            maxThreads = Runtime.getRuntime().availableProcessors();
        }

        System.err.println("Starting kmer mapping at " + new Date());
        System.err.println("*  Number of threads:       " + maxThreads);
        System.err.println("*  Max outstanding tasks:   " + maxTasks);
        System.err.println("*  Nucleotide contigs file: " + nuclContigs);
        System.err.println("*  Reads file:              " + readsFile);
        System.err.println("*  Kmer length:             " + k);

        long startTime = System.currentTimeMillis();
        final ReadExtractor kmerCounter = new ReadExtractor(nuclContigs, k);
        System.err.println("Kmer trie built in " + (System.currentTimeMillis() - startTime) + " ms");

        FastaWriter out = new FastaWriter(System.out);


        SequenceReader reader = new SequenceReader(readsFile);
        Sequence seq;

        startTime = System.currentTimeMillis();
        while ((seq = reader.readNextSequence()) != null) {
            kmerCounter.processRead(seq, out);
        }

        reader.close();

        //service.shutdown();
        //service.awaitTermination(1, TimeUnit.DAYS);

        System.err.println("Processed " + kmerCounter.processedSeqs + " and wrote " + kmerCounter.writtenSeqs + " sequences in " + (System.currentTimeMillis() - startTime) + " ms");

        //System.err.println("Unique kmers in contigs: " + kmerCounter.trie.uniqueWords());
        System.err.println("Processing complete");
    }
}
