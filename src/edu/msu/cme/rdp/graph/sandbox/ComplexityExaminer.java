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
package edu.msu.cme.rdp.graph.sandbox;

import edu.msu.cme.rdp.graph.filter.BloomFilter;
import edu.msu.cme.rdp.graph.filter.CodonWalker;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author fishjord
 */
public class ComplexityExaminer {

    private final BloomFilter bloom;
    private final int maxDepth;
    private final long[] nodesAtDepth;
    private final long[] branchesAtDepth;
    private final long[][] nuclAtDepth;
    private long totalPaths;
    private long totalPathLengths;
    private long totalNodes;
    private long time;

    public ComplexityExaminer(File bloomFile, int maxNuclDepth) throws IOException {
        bloom = BloomFilter.fromFile(bloomFile);
        maxDepth = maxNuclDepth;
        nodesAtDepth = new long[maxDepth];
        branchesAtDepth = new long[maxDepth];
        nuclAtDepth = new long[maxDepth][4];
    }

    public void examine(String kmer) {
        System.out.print(kmer + " left ");
        long startTime = System.currentTimeMillis();
        examine(bloom.new LeftCodonFacade(kmer));
        System.out.println((System.currentTimeMillis() - startTime) / 1000.0 + "s");
        printDetails();

        System.out.print(kmer + " right");
        startTime = System.currentTimeMillis();
        examine(bloom.new RightCodonFacade(kmer));
        System.out.println((System.currentTimeMillis() - startTime) / 1000.0 + "s");
        printDetails();
    }

    private void printDetails() {
        long totalNodes = 0;
        long totalBranches = 0;
        StringBuilder nodeDepth = new StringBuilder();
        StringBuilder branchDepth = new StringBuilder();

        for(int index = 0;index < maxDepth;index++) {
            nodeDepth.append(nodesAtDepth[index]).append("\t");
            branchDepth.append((float)branchesAtDepth[index] / nodesAtDepth[index]).append("\t");

            totalNodes += nodesAtDepth[index];
            totalBranches += branchesAtDepth[index];
        }

        System.out.println("Total nodes:\t" + totalNodes);
        System.out.println("Total Branches:\t" + totalBranches);
        System.out.println("Total paths:\t" + totalPaths);
        System.out.println("Total path length:\t" + totalPathLengths);
        System.out.println();
        System.out.println("Average branching factor:\t" + (double)totalBranches / totalNodes);
        System.out.println("Average path length:\t" + (double)totalPathLengths / totalPaths);
        System.out.println("Nodes at depth:\t" + nodeDepth);
        System.out.println("Branches at depth:\t" + branchDepth);

    }

    private void reset() {
        Arrays.fill(nodesAtDepth, 0);
        Arrays.fill(branchesAtDepth, 0);
        for (int index = 0; index < maxDepth; index++) {
            Arrays.fill(nuclAtDepth[index], 0);
        }
        totalPaths = 0;
        totalPathLengths = 0;
        totalNodes = 0;
        time = System.currentTimeMillis();
    }

    private void examine(CodonWalker walker) {
        reset();
        dfs(walker, 0);
    }

    private void dfs(CodonWalker walker, int depth) {
        Byte b = walker.getNextNucl();
        if (depth >= maxDepth || b == null) {
            totalPaths++;
            totalPathLengths += depth;
            return;
        }

        totalNodes++;
        nodesAtDepth[depth]++;

        if(totalNodes % 1000000 == 0) {
            System.err.println("Processed " + totalNodes + " nodes in " + (System.currentTimeMillis() - time) / 1000.0 + "s");
        }

        do {
            branchesAtDepth[depth]++;
            nuclAtDepth[depth][b]++;

            dfs(walker, depth + 1);
        } while ((b = walker.getSibNucl()) != null);
    }

    public static void main(String[] args) throws IOException {
        if (args.length != 3) {
            System.err.println("USAGE: ComplexityExaminer <bloom_filter> <kmer_file> <max depth>");
            System.exit(1);
        }
        File bloomFile = new File(args[0]);
        Integer maxDepth = Integer.valueOf(args[2]);
        String line;

        final ComplexityExaminer examiner = new ComplexityExaminer(bloomFile, maxDepth);
        Runtime.getRuntime().addShutdownHook(new Thread() {
            public void run() {
                System.err.println("----JAVA SHUTDOWN----");
                examiner.printDetails();
            }
        });

        BufferedReader reader = new BufferedReader(new FileReader(args[1]));
        while((line = reader.readLine()) != null) {
            line = line.trim();
            if(line.isEmpty()) {
                continue;
            }
            examiner.examine(line);
        }
        reader.close();
    }
}
