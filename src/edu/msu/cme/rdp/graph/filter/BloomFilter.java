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
package edu.msu.cme.rdp.graph.filter;

import edu.msu.cme.rdp.graph.hash.CyclicHash;
import edu.msu.cme.rdp.graph.hash.Hash;
import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.readseq.utils.NuclBinMapping;
import java.io.*;
import java.util.BitSet;
import java.util.Date;

/**
 *
 * @author wangqion This is not thread safe
 */
public class BloomFilter implements Serializable {

    private static final int MAX_ASCII = 128;
    private static final int LONGSIZE = 64;
    private static final int MAX_BITSETSIZELOG2 = 30;
    public static final long serialVersionUID = -8788171152437524877L;
    private final Hash hasher;
    private final BitSet[] bitsetArray;
    /**
     * ***********
     * These variables define our bloom filter they can never change once a
     * bloom filter has been created
     ************
     */
    private final int bitsetSizeLog2; // the size of one bitSet 2^bitsetSizeLog2
    private final int hashSizeLog2; // log2 of hash size
    private final long bitsetMask;
    private final long hashMask;
    private final int hashCount;   // number of hash functions
    private final int kmerSize;  // should be less than 32
    private final int bitsetSize;  // the size of one bitSet
    /**
     * ***********
     * These variables keep track of stats about what is in the bloom filter
     ************
     */
    private long uniqueKmers;
    private long totalKmers;
    private long totalStrings;
    private final Date createdOn;

    public static BloomFilter fromFile(File f) throws IOException {
        ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(f)));
        try {
            BloomFilter ret = (BloomFilter) ois.readObject();
            return ret;
        } catch (ClassNotFoundException e) {
            throw new RuntimeException(e);
        } finally {
            ois.close();
        }
    }

    /**
     *
     */
    public BloomFilter(int hashSizeLog2, int hashCount, int kmerSize, int bitsetSizeLog2) {
        createdOn = new Date();
        this.hashCount = hashCount;
        this.hashSizeLog2 = hashSizeLog2;
        this.hashMask = (1L << hashSizeLog2) - 1;
        this.kmerSize = kmerSize;
        this.bitsetSizeLog2 = bitsetSizeLog2;

        if (bitsetSizeLog2 > MAX_BITSETSIZELOG2) {
            throw new IllegalArgumentException("Can't have a bitset larger than 2^" + MAX_BITSETSIZELOG2);
        }
        if (hashSizeLog2 > LONGSIZE) {
            throw new IllegalArgumentException("Can't have filter larger than 2^" + LONGSIZE);
        }
        this.bitsetSize = (1 << this.bitsetSizeLog2);
        this.bitsetMask = this.bitsetSize - 1;
        int bitSetCount = 1;

        if (hashSizeLog2 > bitsetSizeLog2) {
            bitSetCount = (1 << (hashSizeLog2 - bitsetSizeLog2));
        }


        hasher = new CyclicHash(kmerSize);
        bitsetArray = new BitSet[bitSetCount];
        for (int i = 0; i < bitSetCount; i++) {
            bitsetArray[i] = new BitSet(bitsetSize);
        }

    }


    /*
     * Returns True if add was successful, and kmer was not already present Sets
     * internal state to point to this node
     */
    boolean addNode(long xHash, long yHash) {
        boolean wasSet = false;
        //to avoid overflow, we get the mod first because (a+B)%m = (a%m + b%m )%m
        //Only works if hash size < maxLong in size.
        xHash = xHash & hashMask;
        yHash = yHash & hashMask;
        long hashVal = xHash;
        for (int i = 0; i < hashCount; ++i) {
            wasSet |= setBit(hashVal);
            hashVal += yHash;
            hashVal &= hashMask;
        }
        return wasSet;
    }

    boolean hasNode(long xHash, long yHash) {
        boolean wasSet = true;
        //to avoid overflow, we get the mod first because (a+B)%m = (a%m + b%m )%m
        //Only works if hash size < maxLong in size.
        xHash = xHash & hashMask;
        yHash = yHash & hashMask;
        long hashVal = xHash;
        for (int i = 0; i < hashCount && wasSet; ++i) {
            wasSet &= isSet(hashVal);
            hashVal += yHash;
            hashVal &= hashMask;
        }
        return wasSet;
    }

    /*
     * Returns True if bit not previously set
     */
    boolean setBit(long bit) {
        //bit = bit & hashMask;
        int bitsetOffset = (int) (bit & bitsetMask);
        int bitsetNo = (int) (bit >>> bitsetSizeLog2);

        boolean wasSet = bitsetArray[bitsetNo].get(bitsetOffset);
        bitsetArray[bitsetNo].set(bitsetOffset);
        return !wasSet;
    }

    boolean isSet(long bit) {
        bit = bit & hashMask;
        int bitsetOffset = (int) (bit & bitsetMask);
        int bitsetNo = (int) (bit >> bitsetSizeLog2);
        return bitsetArray[bitsetNo].get(bitsetOffset);
    }

    public int getBitsetSize() {
        return bitsetSize;
    }

    public int getHashCount() {
        return hashCount;
    }

    public int getHashSizeLog2() {
        return hashSizeLog2;
    }

    public int getKmerSize() {
        return kmerSize;
    }

    public long getBitsetMask() {
        return bitsetMask;
    }

    public int getBitsetSizeLog2() {
        return bitsetSizeLog2;
    }

    public long getHashMask() {
        return hashMask;
    }

    public String getHasherClassName() {
        return hasher.getClass().getCanonicalName();
    }

    public long getTotalKmers() {
        return totalKmers;
    }

    public long getTotalStrings() {
        return totalStrings;
    }

    public long getUniqueKmers() {
        return uniqueKmers;
    }

    public Date getCreatedOn() {
        return createdOn;
    }

    public int getNumBitsets() {
        return bitsetArray.length;
    }
    public static byte[] next = new byte[4];

    static {
        next[NuclBinMapping.a] = NuclBinMapping.t;
        next[NuclBinMapping.t] = NuclBinMapping.g;
        next[NuclBinMapping.g] = NuclBinMapping.c;
        next[NuclBinMapping.c] = NuclBinMapping.a;
    }

    public class GraphState {

        protected long rcHashValue;
        protected long fwdHashValue;
        protected byte[] kmer = new byte[kmerSize];
        protected byte[] rkmer = new byte[kmerSize];
        protected int kmerLeftIdx, rkmerLeftIdx;

        public GraphState(char[] s) {
            setState(s);
        }

        public GraphState() {
        }

        public void setState(char[] s) {
            if (s.length != kmerSize) {
                throw new InvalidDNABaseException("input length not equal to k-mer length");
            }

            for (int i = 0; i < kmerSize; ++i) {
                byte c = NuclBinMapping.validateLookup[s[i]];
                if (c == -1) {
                    throw new InvalidDNABaseException("Input contains non nucleotide character: " + c);
                }
                kmer[i] = c;
                rkmer[kmerSize - 1 - i] = NuclBinMapping.complementLookup[c];
            }
            fwdHashValue = 0;
            rcHashValue = 0;
            for (int i = 0; i < kmerSize; ++i) {
                fwdHashValue = hasher.eatRight(fwdHashValue, kmer[i]);
                rcHashValue = hasher.eatRight(rcHashValue, rkmer[i]);
            }
            kmerLeftIdx = 0;
            rkmerLeftIdx = 0;
        }
        int fwd;
        int rc;

        public void setStateDangerously(byte[] s) {
            fwdHashValue = 0;
            rcHashValue = 0;
            for (int i = 0; i < kmerSize; ++i) {
                fwd = s[i];
                rc = ((fwd == NuclBinMapping.a)? NuclBinMapping.t : ((fwd == NuclBinMapping.c)? NuclBinMapping.g : ((fwd == NuclBinMapping.g)? NuclBinMapping.c : NuclBinMapping.a)));
                //kmer[i] = s[i];
                //rkmer[kmerSize - 1 - i] = rc;
                fwdHashValue = hasher.eatRight(fwdHashValue, fwd);
                rcHashValue = hasher.eatRight(rcHashValue, rc);
            }
            kmerLeftIdx = 0;
            rkmerLeftIdx = 0;
        }

        //protected void clearState() {
        public void clearState() {
            fwdHashValue = rcHashValue = 0;
            kmerLeftIdx = rkmerLeftIdx = 0;
            for (int i = 0; i < kmerSize; ++i) {
                kmer[i] = rkmer[i] = 0;
            }

        }

        public boolean hasCurrent() {
            long xHash = (fwdHashValue > rcHashValue) ? fwdHashValue : rcHashValue;
            long yHash = (fwdHashValue > rcHashValue) ? rcHashValue : fwdHashValue;

            return BloomFilter.this.hasNode(xHash, yHash);
        }

        //protected void loadCharRight( char inChar) {
        public void loadCharRight(char inChar) {
            byte c = NuclBinMapping.validateLookup[inChar];
            if (c == -1) {
                throw new InvalidDNABaseException("Input contains non nucleotide character: " + inChar);
            }
            fwdHashValue = hasher.eatRight(fwdHashValue, c);
            kmer[kmerLeftIdx] = c;
            if (++kmerLeftIdx >= kmerSize) {
                kmerLeftIdx -= kmerSize;
            }

            c = NuclBinMapping.complementLookup[c];
            if (--rkmerLeftIdx < 0) {
                rkmerLeftIdx += kmerSize;
            }
            rcHashValue = hasher.eatLeft(rcHashValue, c);
            rkmer[rkmerLeftIdx] = c;
        }

        /**
         * append a char to the left, update the state
         *
         * @param c
         */
        //protected void shiftRight( char inChar) {
        public void shiftRight(char inChar) {
            byte c = NuclBinMapping.validateLookup[inChar];
            if (c == -1) {
                throw new InvalidDNABaseException("Input contains non nucleotide character: " + inChar);
            }
            fwdHashValue = hasher.updateRight(fwdHashValue, kmer[kmerLeftIdx], c);
            kmer[kmerLeftIdx] = c;
            if (++kmerLeftIdx >= kmerSize) {
                kmerLeftIdx -= kmerSize;
            }

            c = NuclBinMapping.complementLookup[c];
            if (--rkmerLeftIdx < 0) {
                rkmerLeftIdx += kmerSize;
            }
            rcHashValue = hasher.updateLeft(rcHashValue, rkmer[rkmerLeftIdx], c);
            rkmer[rkmerLeftIdx] = c;
        }

        /**
         * append a char to the right, update the state
         */
        protected void shiftLeft(char inChar) {
            byte c = NuclBinMapping.validateLookup[inChar];
            if (c == -1) {
                throw new InvalidDNABaseException("Input contains non nucleotide character: " + inChar);
            }
            if (--kmerLeftIdx < 0) {
                kmerLeftIdx += kmerSize;
            }
            fwdHashValue = hasher.updateLeft(fwdHashValue, kmer[kmerLeftIdx], c);
            kmer[kmerLeftIdx] = c;

            c = NuclBinMapping.complementLookup[c];
            rcHashValue = hasher.updateRight(rcHashValue, rkmer[rkmerLeftIdx], c);
            rkmer[rkmerLeftIdx] = c;
            if (++rkmerLeftIdx >= kmerSize) {
                rkmerLeftIdx -= kmerSize;
            }
        }
    }

    public class GraphBuilder extends BloomFilter.GraphState {

        private long numStr = 0;
        private long numkmer = 0;
        private long numUniqueKmer = 0;

        public GraphBuilder() {
        }

        /**
         * @return true if bloomfilter has the kmer represented by the current
         * hashvalues
         */
        public boolean setCurrent() {
            long xHash = (fwdHashValue > rcHashValue) ? fwdHashValue : rcHashValue;
            long yHash = (fwdHashValue > rcHashValue) ? rcHashValue : fwdHashValue;
            return BloomFilter.this.addNode(xHash, yHash);
        }

        /**
         * adds all the kmers from the input seqStr to bloomfilter will ignore
         * the kmers with invalid bases
         *
         * @param seqStr
         */
        public void addString(char[] seqStr) {
            numStr++;
            BloomFilter.this.totalStrings++;

            int i = 0;
            while (i < seqStr.length) {
                try {
                    clearState();
                    int j;
                    for (j = 0; j < kmerSize && i < seqStr.length; ++j) {
                        loadCharRight(seqStr[i]);
                        ++i;
                    }
                    if (j < kmerSize) {
                        break;
                    }


                    boolean wasSet = setCurrent();
                    numkmer++;
                    numUniqueKmer += wasSet ? 1 : 0;

                    while (i < seqStr.length) {
                        shiftRight(seqStr[i]);
                        ++i;
                        wasSet = setCurrent();
                        numkmer++;
                        BloomFilter.this.totalKmers++;

                        numUniqueKmer += wasSet ? 1 : 0;
                        BloomFilter.this.uniqueKmers += wasSet ? 1 : 0;
                    }

                } catch (InvalidDNABaseException e) {
                    ++i;
                }
            } // end while
        }

        /**
         *
         * @return number of kmers attempted to added to bloomfilter
         */
        public long getKmerAdded() {
            return numkmer;
        }

        /**
         *
         * @return nubmer of unique kmers added to bloomfilter
         */
        public long getUniqueKmerAdded() {
            return numUniqueKmer;
        }
    }

    /**
     * a class to get codons from the right side of the starting kmer
     */
    public abstract class CodonFacade implements CodonWalker {

        public long rcHashValue = 0;
        public long fwdHashValue = 0;
        //protected StringBuilder path = new StringBuilder();
        public int framePtr = -1;      // index of the end of last codon
        private static final byte startChar = NuclBinMapping.a;
        /*
         * The PathHolder imposes a constraint on how the path is represented
         * The 'left most' character in a path (ie path[0]) is always the earliest
         * one pushed on, and the right most is the most recent pushed on.
         *
         * This is used to simplify bitwise operations (shift the kmer, | the new emission)
         */
        protected PathHolder path = new PathHolder();
        protected int pathPtr = -1;

        public CodonFacade(String s) {
            this(s.toCharArray());
        }

        public CodonFacade(char[] s) {
            jumpTo(s);
        }

        public long getFwdHash() {
            return fwdHashValue;
        }

        public long getRcHash() {
            return rcHashValue;
        }

        public final void jumpTo(char[] kmer) {
            if (kmer.length != kmerSize) {
                throw new InvalidDNABaseException("input length [" + kmer.length + "] not equal to k-mer length[" + kmerSize + "]: " + new String(kmer));
            }

            initialize(kmer);

            if (!hasCurrent()) {
                throw new InvalidDNABaseException("kmer not in bloomfilter: " + new String(kmer));
            }
        }

        public final void jumpTo(Kmer kmer, long fwdHash, long rcHash) {
            reset(kmer, fwdHash, rcHash);

            if (!hasCurrent()) {
                throw new IllegalArgumentException("kmer not in bloomfilter: " + kmer);
            }
        }

        protected abstract void initialize(char[] s);

        /**
         * does not change state.
         *
         * @return true if the bloomfilter has the kmer represented by the
         * current hashvalues
         */
        protected final boolean hasCurrent() {
            long xHash = (fwdHashValue > rcHashValue) ? fwdHashValue : rcHashValue;
            long yHash = (fwdHashValue > rcHashValue) ? rcHashValue : fwdHashValue;

            return BloomFilter.this.hasNode(xHash, yHash);
        }

        protected final void reset(Kmer kmer, long fwdHash, long rcHash) {
            fwdHashValue = fwdHash;
            rcHashValue = rcHash;
            //initialize(kmer);
            framePtr = kmerSize - 1;// - (kmerSize % 3);
            pathPtr = kmerSize - 1;
            path.init(kmer);
        }

        protected abstract void updateHashForward(byte out, byte in);

        protected abstract void updateHashReverse(byte out, byte in);

        protected abstract void replaceHash(byte out, byte in);

        /**
         * Returns amino acid for next codon, if one exists. Updates state to
         * new codon. If none exists, does not change state.
         *
         * @return
         */
        public abstract NextCodon getNextCodon();

        /**
         * Recursive call. Loops over all existing vertices to the right and
         * calls self. Returns true when a full codon is found with state set to
         * new codon. returns false and does not change state if new codon not
         * found.
         *
         * @return
         */
        protected boolean finishCodon() {
            if ((pathPtr - framePtr) == 3) {
                return true;
            }
            if (probe()) {  // increments pathPtr
                do {
                    if (finishCodon()) {
                        return true;
                    }
                } while (replace()); // returns false if choices exausted
                //couldn't complete codon
                backup();
            }
            return false;
        }

        protected boolean walk(byte... walk) {
            for (byte w : walk) {
                pathPtr++;
                path.push(w);

                byte leftChar = path.get(pathPtr - kmerSize);
                updateHashForward(leftChar, w);
            }

            if (hasCurrent()) {
                return true;
            } else {
                return false;
            }
        }

        protected boolean walkUpdate(byte... walk) {
            boolean ret = walk(walk);
            framePtr = pathPtr;

            return ret;
        }

        /**
         * @return alternate amino acid for current codon position. If fails,
         * removes current codon (previous codon becomes current). If it backs
         * up to the starting k-mer, it leaves it in that state and returns 0
         *
         */
        public NextCodon getSibCodon() {
            NextCodon retVal = null;
            framePtr = framePtr - 3;
            while (framePtr < pathPtr && pathPtr >= kmerSize) {
                while (retVal == null && replace()) {
                    retVal = getNextCodon();
                    if (retVal != null) {
                        return retVal;
                    }
                }
                backup();
            }

            return retVal;
        }

        public boolean hasMoreCodons() {
            return (pathPtr >= kmerSize);
        }

        /**
         * attempts to find a vertex with right-most character replaced, in
         * order defined in next[] if found, set state to that vertex and
         * returns true, repeated calls to this method will find all such
         * vertices. if not found, returns false and does not change state
         *
         * @return
         */
        protected final boolean replace() {
            byte origChar = path.get(pathPtr);
            byte oldC = origChar;
            byte newC = next[oldC];

            while (newC != 0) {
                replaceHash(oldC, newC);
                if (hasCurrent()) {
                    path.remove();
                    path.push(newC);
                    return true;
                }
                oldC = newC;
                newC = next[oldC];
            }

            replaceHash(oldC, origChar);
            // didn't find a good replacement character, restore to the starting state

            return false;
        }

        /**
         * move pathptr one char left, update the state
         */
        private void backup() {
            byte newC = path.get(pathPtr - kmerSize);
            byte oldC = path.remove();
            pathPtr--;

            updateHashReverse(oldC, newC);
        }

        /**
         *
         * if a vertex to the right exists, move the pathPtr to one char right,
         * sets internal state to that vertex and return true. if a vertex to
         * the right does not exist, do not change state and return false;
         */
        protected boolean probe() {
            if (walk(startChar)) {
                return true;
            } else if (replace()) {
                return true;
            } else {
                backup();
                return false;
            }
        }

        /**
         *
         * @return the path starting from the char right after the original
         * kmer, all the way to the current char pointed by pathPtr
         */
        public String getPathString() {
            return new String(path.toCharArray()).substring(kmerSize);
        }

        public Byte getNextNucl() {
            if (!probe()) {
                return null;
            }

            return path.get(pathPtr);
        }

        public Byte getSibNucl() {
            if (replace()) {
                return path.get(pathPtr);
            }

            backup();
            return null;
        }

        public boolean hasMoreNucl() {
            return hasMoreCodons();
        }

	public int getLength() {
	    return pathPtr - kmerSize;
	}
    }

    /**
     * a class to get codons from the right side of the starting kmer
     */
    public final class RightCodonFacade extends CodonFacade {

        public RightCodonFacade(String s) {
            super(s);
        }

        public RightCodonFacade(char[] s) {
            super(s);
        }

        protected void initialize(char[] s) {
            fwdHashValue = 0;
            rcHashValue = 0;

            for (int i = 0; i < kmerSize; ++i) {
                byte c = NuclBinMapping.validateLookup[s[i]];
                if (c == -1) {
                    throw new InvalidDNABaseException("Input contains non nucleotide character: " + s[i]);
                }
                fwdHashValue = hasher.eatRight(fwdHashValue, c);
                rcHashValue = hasher.eatLeft(rcHashValue, NuclBinMapping.complementLookup[c]);
            }

            path.init(new Kmer(s));
            framePtr = kmerSize - 1 - (kmerSize % 3);
            pathPtr = kmerSize - 1;
        }

        protected void replaceHash(byte out, byte in) {
            fwdHashValue = hasher.replaceRight(fwdHashValue, out, in);
            rcHashValue = hasher.replaceLeft(rcHashValue, NuclBinMapping.complementLookup[out], NuclBinMapping.complementLookup[in]);
        }

        protected void updateHashForward(byte out, byte in) {
            fwdHashValue = hasher.updateRight(fwdHashValue, out, in);
            rcHashValue = hasher.updateLeft(rcHashValue, NuclBinMapping.complementLookup[out], NuclBinMapping.complementLookup[in]);
        }

        protected void updateHashReverse(byte out, byte in) {
            fwdHashValue = hasher.updateLeft(fwdHashValue, out, in);
            rcHashValue = hasher.updateRight(rcHashValue, NuclBinMapping.complementLookup[out], NuclBinMapping.complementLookup[in]);
        }

        public NextCodon getNextCodon() {
            if (!finishCodon()) {
                return null;
            }
            framePtr = pathPtr;

            int i = path.size();
            return new NextCodon(true, path.get(i - 3), path.get(i - 2), path.get(i - 1));
        }
    }

    /**
     * a class to get codons from the right side of the starting kmer
     * To comply with the constraints set by the PathHolder we store the
     * left search path in reverse order.  This has several implications, including
     * that codons are stored in reverse order
     */
    public final class LeftCodonFacade extends CodonFacade {

        public LeftCodonFacade(String s) {
            super(s);
        }

        public LeftCodonFacade(char[] s) {
            super(s);
        }

        protected void initialize(char[] s) {
            fwdHashValue = 0;
            rcHashValue = 0;

            /*
             * Be cause we're storing the path in reverse order we need to
             * load the initial kmer in the reverse order
             */
            for (int i = kmerSize - 1; i >= 0; i--) {
                byte c = NuclBinMapping.validateLookup[s[i]];
                if (c == -1) {
                    throw new InvalidDNABaseException("Input contains non nucleotide character: " + s[i]);
                }
                fwdHashValue = hasher.eatLeft(fwdHashValue, c);
                rcHashValue = hasher.eatRight(rcHashValue, NuclBinMapping.complementLookup[c]);
            }

            path.init(new Kmer(new StringBuilder(new String(s)).reverse().toString().toCharArray()));
            framePtr = kmerSize - 1 - (kmerSize % 3);
            pathPtr = kmerSize - 1;
        }

        protected void replaceHash(byte out, byte in) {
            fwdHashValue = hasher.replaceLeft(fwdHashValue, out, in);
            rcHashValue = hasher.replaceRight(rcHashValue, NuclBinMapping.complementLookup[out], NuclBinMapping.complementLookup[in]);
        }

        protected void updateHashForward(byte out, byte in) {
            fwdHashValue = hasher.updateLeft(fwdHashValue, out, in);
            rcHashValue = hasher.updateRight(rcHashValue, NuclBinMapping.complementLookup[out], NuclBinMapping.complementLookup[in]);
        }

        protected void updateHashReverse(byte out, byte in) {
            fwdHashValue = hasher.updateRight(fwdHashValue, out, in);
            rcHashValue = hasher.updateLeft(rcHashValue, NuclBinMapping.complementLookup[out], NuclBinMapping.complementLookup[in]);
        }

        public NextCodon getNextCodon() {
            if (!finishCodon()) {
                return null;
            }
            framePtr = pathPtr;

            int i = path.size();
            //Yes, we're using the last three from the kmer, an it does look
            //like this is the wrong order, but NextCodon needs them in the
            //'Graph walk order', the first argument tells the constructor
            //that it needs to reverse the codon when translating to protein
            return new NextCodon(false, path.get(i - 3), path.get(i - 2), path.get(i - 1));
        }

        @Override
        public String getPathString() {
            //Since we're storing the left side path in reverese order, reverse
            //it before we return it
            return new StringBuilder(super.getPathString()).reverse().toString();
        }
    }
}
