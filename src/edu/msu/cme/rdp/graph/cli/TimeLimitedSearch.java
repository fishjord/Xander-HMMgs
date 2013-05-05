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
import edu.msu.cme.rdp.graph.search.HMMGraphSearch;
import edu.msu.cme.rdp.graph.search.HMMGraphSearch.HackTerminateException;
import edu.msu.cme.rdp.graph.search.SearchResult;
import edu.msu.cme.rdp.graph.search.SearchTarget;
import edu.msu.cme.rdp.kmer.io.KmerStart;
import edu.msu.cme.rdp.kmer.io.KmerStartsReader;
import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.FutureTask;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

/**
 *
 * @author fishjord
 */
public class TimeLimitedSearch {

    private static class TimeLimitedSearchThread implements Callable<List<SearchResult>> {

        private HMMGraphSearch searchMethod;
        private SearchTarget target;

        public TimeLimitedSearchThread(HMMGraphSearch searchMethod, SearchTarget target) {
            this.searchMethod = searchMethod;
            this.target = target;
        }

        public List<SearchResult> call() throws Exception {
            try {
                return searchMethod.search(target);
            } catch (HackTerminateException e) {
                return null;
            }
        }
    }

    private static class TimeStamppedFutureTask extends FutureTask<List<SearchResult>> {

        private long startedAt = -1;
        private String startingWord;

        public TimeStamppedFutureTask(Runnable runnable, List<SearchResult> result) {
            super(runnable, result);
        }

        public TimeStamppedFutureTask(Callable callable) {
            super(callable);
        }

        public TimeStamppedFutureTask(TimeLimitedSearchThread callable) {
            super(callable);
            startingWord = callable.target.getKmer();
        }

        @Override
        public void run() {
            this.startedAt = System.currentTimeMillis();
            super.run();
        }

        public boolean hasStarted() {
            return startedAt != -1;
        }

        public long getStartedAt() {
            return startedAt;
        }

        public String getStartingWord() {
            return startingWord;
        }

        @Override
        public boolean cancel(boolean c) {
            return super.cancel(c);
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 6 && !(args.length == 7 && args[0].equals("-u"))) {
            System.err.println("USAGE: TimeLimitedSearch -u <k> <limit_in_seconds> <bloom_filter> <for_hmm> <rev_hmm> <kmers>");
            System.exit(1);
        }

        boolean normalized = true;
        if(args.length == 7) {
            normalized = false;
            args = Arrays.copyOfRange(args, 1, args.length);
        }

        int k = Integer.valueOf(args[0]);
        long timeLimit = Long.valueOf(args[1]);

        File bloomFile = new File(args[2]);
        File forHMMFile = new File(args[3]);
        File revHMMFile = new File(args[4]);
        File kmersFile = new File(args[5]);

        File nuclOutFile = new File(kmersFile.getName() + "_nucl.fasta");
        File alignOutFile = new File(kmersFile.getName() + ".alignment");
        File protOutFile = new File(kmersFile.getName() + "_prot.fasta");

        HMMGraphSearch search = new HMMGraphSearch(k);

        ProfileHMM forHMM;
        ProfileHMM revHMM;

        if(normalized) {
             forHMM = HMMER3bParser.readModel(forHMMFile);
             revHMM = HMMER3bParser.readModel(revHMMFile);
        } else {
             forHMM = HMMER3bParser.readUnnormalized(forHMMFile);
             revHMM = HMMER3bParser.readUnnormalized(revHMMFile);
        }

        FastaWriter nuclOut = new FastaWriter(nuclOutFile);
        FastaWriter alignOut = new FastaWriter(alignOutFile);
        FastaWriter protOut = null;
        boolean isProt = forHMM.getAlphabet() == SequenceType.Protein;

        if (isProt) {
            protOut = new FastaWriter(protOutFile);
        }

        int kmerCount = 0;
        int contigCount = 1;

        long startTime;

        startTime = System.currentTimeMillis();
        ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(bloomFile)));
        BloomFilter bloom = (BloomFilter) ois.readObject();
        ois.close();
        System.err.println("Bloom filter loaded in " + (System.currentTimeMillis() - startTime) + " ms");

        System.err.println("Starting hmmgs search at " + new Date());
        System.err.println("*  Kmer file:               " + kmersFile);
        System.err.println("*  Bloom file:              " + bloomFile);
        System.err.println("*  Forward hmm file:        " + forHMMFile);
        System.err.println("*  Reverse hmm file:        " + revHMMFile);
        System.err.println("*  Searching prot?:         " + isProt);
        System.err.println("*  # paths:                 " + k);
        System.err.println("*  Nucl contigs out file    " + nuclOutFile);
        System.err.println("*  Prot contigs out file    " + protOutFile);

        startTime = System.currentTimeMillis();
        HMMBloomSearch.printHeader(System.out, isProt);

        //Set<String> processed = new HashSet();
        String key;

        KmerStart line;
        KmerStartsReader reader = new KmerStartsReader(kmersFile);

        try {
            while ((line = reader.readNext()) != null) {

                kmerCount++;

                if (line.getMpos() == 0) {
                    System.err.println("Skipping line " + line);
                    continue;
                }

                TimeStamppedFutureTask future = new TimeStamppedFutureTask(
                        new TimeLimitedSearchThread(search,
                        new SearchTarget(line.getGeneName(),
                        line.getQueryId(), line.getRefId(), line.getNuclKmer(), 0,
                        line.getMpos() - 1, forHMM, revHMM, bloom)));

                Thread t = new Thread(future);
                t.setDaemon(true);
                t.setPriority(Thread.MAX_PRIORITY);
                t.start();

                try {
                    List<SearchResult> searchResults = future.get(timeLimit, TimeUnit.SECONDS);

                    for (SearchResult result : searchResults) {
                        String seqid = "contig_" + (contigCount++);

                        HMMBloomSearch.printResult(seqid, isProt, result, System.out);

                        nuclOut.writeSeq(seqid, result.getNuclSeq());
                        alignOut.writeSeq(seqid, result.getAlignSeq());
                        if (isProt) {
                            protOut.writeSeq(seqid, result.getProtSeq());
                        }
                    }

                } catch (TimeoutException e) {
                    System.out.println("-\t" + future.getStartingWord() + (isProt ? "\t-" : "") + "\t-\t-\t-\t-");
                    future.cancel(true);
                } catch (Exception e) {
                    System.out.println("-\t" + future.getStartingWord() + (isProt ? "\t-" : "") + "\t-\t-\t-\t-");
                    e.printStackTrace();
                    if (e.getCause() != null) {
                        e.getCause().printStackTrace();
                    }
                    future.cancel(true);
                }
            }
            System.err.println("Read in " + kmerCount + " kmers and created " + contigCount + " contigs in " + (System.currentTimeMillis() - startTime) / 1000f + " seconds");
        } finally {
            nuclOut.close();
            if (isProt) {
                protOut.close();
            }
            System.out.close();
        }

    }
}
