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

import edu.msu.cme.rdp.graph.utils.BloomFilterAppender;
import edu.msu.cme.rdp.graph.utils.ContigMerger;
import edu.msu.cme.rdp.graph.utils.BloomFilterStats;
import java.util.Arrays;

/**
 *
 * @author fishjord
 */
public class HMMgsMain {

    private static void printLicense() {
        System.out.println("Xander  Copyright (C) 2012  Michigan State University\n"
                + "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.\n"
                + "This is free software, and you are welcome to redistribute it\n"
                + "under certain conditions; see LICENSE for details");
    }

    private static void printUsageAndExit() {
        System.err.println("USAGE: HMMgs <command> <options>");
        System.err.println("\tbuild       - Build a bloom filter");
        System.err.println("\tstats       - Display bloom filter stats");
        System.err.println("\tsearch      - Search a bloom filter with an hmm");
        System.err.println("\tmerge       - Merge HMMgs left and right fragments");
        System.err.println("\tlicense     - Print the license");
        System.exit(1);
    }

    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            printUsageAndExit();
        }

        String cmd = args[0];
        args = Arrays.copyOfRange(args, 1, args.length);

        if (cmd.equals("build")) {
            BloomFilterBuilder.main(args);
        } else if (cmd.equals("append")) {
            BloomFilterAppender.main(args);
        } else if (cmd.equals("stats")) {
            BloomFilterStats.main(args);
        } else if (cmd.equals("search")) {
            TimeLimitedSearch.main(args);
        } else if (cmd.equals("basic")) {
            BasicSearch.main(args);
        } else if (cmd.equals("merge")) {
            ContigMerger.main(args);
        } else if (cmd.equals("license")) {
            printLicense();
        } else {
            printUsageAndExit();
        }
    }
}
