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
package edu.msu.cme.rdp.graph.visual;

import edu.msu.cme.rdp.alignment.hmm.HMMER3bParser;
import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import edu.msu.cme.rdp.graph.cli.HMMBloomSearch;
import edu.msu.cme.rdp.graph.filter.BloomFilter;
import edu.msu.cme.rdp.graph.search.AStarNode;
import edu.msu.cme.rdp.graph.search.HMMGraphSearch;
import edu.msu.cme.rdp.graph.search.HMMGraphSearch.HackTerminateException;
import edu.msu.cme.rdp.graph.search.SearchTarget;
import edu.msu.cme.rdp.kmer.io.KmerStart;
import edu.msu.cme.rdp.kmer.io.KmerStartsReader;
import edu.msu.cme.rdp.readseq.SequenceType;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

/**
 *
 * @author fishjord
 */
public class GraphSearch {

    private static class TimeLimitedSearchThread implements Callable<List<AStarNode>> {

        private HMMGraphSearch searchMethod;
        private SearchTarget target;

        public TimeLimitedSearchThread(HMMGraphSearch searchMethod, SearchTarget target) {
            this.searchMethod = searchMethod;
            this.target = target;
        }

        public List<AStarNode> call() throws Exception {
            try {
                return searchMethod.searchGraph(target);
            } catch (HackTerminateException e) {
                return null;
            }
        }
    }

    private static class TimeStamppedFutureTask extends FutureTask<List<AStarNode>> {

        private long startedAt = -1;
        private String startingWord;

        public TimeStamppedFutureTask(Runnable runnable, List<AStarNode> result) {
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
        if (args.length != 6 && args.length != 7) {
            System.err.println("USAGE: TimeLimitedSearch <k> <limit_in_seconds> <bloom_filter> <for_hmm> <rev_hmm> <kmers> [threads=#processors]");
            System.exit(1);
        }

        int k = Integer.valueOf(args[0]);
        long timeLimit = Long.valueOf(args[1]) * 1000;

        File bloomFile = new File(args[2]);
        File forHMMFile = new File(args[3]);
        File revHMMFile = new File(args[4]);
        File kmersFile = new File(args[5]);

        File graphOutStem = new File(kmersFile.getName() + "_graph_");

        HMMGraphSearch search = new HMMGraphSearch(k);

        ProfileHMM forHMM = HMMER3bParser.readModel(forHMMFile);
        ProfileHMM revHMM = HMMER3bParser.readModel(revHMMFile);
        boolean isProt = forHMM.getAlphabet() == SequenceType.Protein;

        int threads = Runtime.getRuntime().availableProcessors();

        if (args.length == 7) {
            threads = Integer.valueOf(args[6]);
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
        System.err.println("*  Number of threads:       " + threads);
        System.err.println("*  Kmer file:               " + kmersFile);
        System.err.println("*  Bloom file:              " + bloomFile);
        System.err.println("*  Forward hmm file:        " + forHMMFile);
        System.err.println("*  Reverse hmm file:        " + revHMMFile);
        System.err.println("*  Searching prot?:         " + isProt);
        System.err.println("*  # paths:                 " + k);
        System.err.println("*  Graph out stem:          " + graphOutStem);


        ExecutorService executor = Executors.newFixedThreadPool(threads);

        startTime = System.currentTimeMillis();
        HMMBloomSearch.printHeader(System.out, isProt);

        List<TimeStamppedFutureTask> tasks = new ArrayList();

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
                        line.getQueryId(), line.getRefId(), line.getKmer(), 0,
                        line.getMpos(), forHMM, revHMM, bloom)));

                executor.execute(future);
                tasks.add(future);
            }

            Graph graph = new Graph();
            int searches = 0;
            for (TimeStamppedFutureTask future : tasks) {
                try {

                    long startWaiting = System.currentTimeMillis();
                    while (!future.hasStarted()) {
                        if (System.currentTimeMillis() - startWaiting > timeLimit) {
                            throw new TimeoutException();
                        }
                    }
                    long delta = timeLimit - (System.currentTimeMillis() - future.getStartedAt());

                    if (delta < 0) {
                        throw new TimeoutException();
                    }

                    List<AStarNode> searchResults = future.get(delta, TimeUnit.MILLISECONDS);

                    System.out.println("Search " + ++searches + " / " + tasks.size() + " done");
                    for (AStarNode result : searchResults) {
                        graph.connectAll(result);
                    }
                } catch (TimeoutException e) {
                    System.out.println("Search " + ++searches + " / " + tasks.size() + " canceled");
                    future.cancel(true);
                } catch (Exception e) {
                    System.out.println("Search " + ++searches + " / " + tasks.size() + " failed: " + e.getMessage());
                    future.cancel(true);
                }
            }

            System.out.println("Graph searching done, writing graph");
            graph.writeDot(graphOutStem);

	    ObjectOutputStream ps = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(kmersFile.getName() + "_graph.ser")));
	    ps.writeObject(graph);
	    ps.close();

            executor.shutdown();
            System.err.println("Awaiting thread temination");
            executor.awaitTermination(1, TimeUnit.DAYS);

            System.err.println("Read in " + kmerCount + " kmers and created " + contigCount + " contigs in " + (System.currentTimeMillis() - startTime) / 1000f + " seconds");
        } finally {
            System.out.close();
        }

    }
}
