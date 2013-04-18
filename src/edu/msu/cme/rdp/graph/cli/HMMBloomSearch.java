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

import edu.msu.cme.rdp.graph.search.SearchResult;
import java.io.PrintStream;
import java.text.DecimalFormat;

/**
 *
 * @author fishjord
 */
public class HMMBloomSearch {

    private static final DecimalFormat format = new DecimalFormat("#.##");
    private static final DecimalFormat threeDig = new DecimalFormat("#.###");

    public static void printHeader(PrintStream out) {
        printHeader(out, false);
    }

    public static void printHeader(PrintStream out, boolean protSearch) {
        out.println("#contig_id\tgene_name\tquery_id\trefseq_id\tstarting kmer\tstarting state\tnucl length" + (protSearch ? "\tprot length" : "") + "\tsearch direction\tnats\tbits\ttime (s)");
    }

    public static void printResult(String seqid, SearchResult result, PrintStream out) {
        printResult(seqid, true, result, out);
    }

    public static void printResult(String seqid, boolean protSearch, SearchResult result, PrintStream out) {
        System.out.println(
                seqid + "\t"
                + result.getStart().getGeneName() + "\t"
                + result.getStart().getQuerySeqid() + "\t"
                + result.getStart().getRefSeqid() + "\t"
                + result.getKmer() + "\t"
		+ result.getModelPosition() + "\t"
                + result.getNuclSeq().length() + "\t"
                + (protSearch ? result.getProtSeq().length() + "\t" : "")
                + result.getSearchDirection() + "\t"
                + ((Double.isInfinite(result.getNats()) || Double.isNaN(result.getNats()))? "?" : threeDig.format(result.getNats())) + "\t"
                + format.format(result.getBits()) + "\t"
                + threeDig.format(result.getTime()  / 1000.0f));
    }
}
