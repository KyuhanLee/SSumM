package ssumm;

import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ThreadLocalRandom;

public class MinH implements Callable<int[]> {
    public int numSuperNodes, numNodes;
    public IntArrayList superList;
    public IntArrayList[] insideSupernode;
    public long[] edges;


    public MinH(int numSuperNodes, int numNodes, IntArrayList superList, IntArrayList[] insideSupernode, long[] edges){
        this.numSuperNodes = numSuperNodes; this.numNodes = numNodes;
        this.superList = superList;
        this.insideSupernode = insideSupernode;
        this.edges = edges;
    }


    public int[] getRandomPermutation (int length){
        int[] array = new int[length];
        ThreadLocalRandom random = ThreadLocalRandom.current();
        for(int i=0;i<length;i++){
            array[i] = i;
        }
        for(int i = 0; i < length; i++){
            int ran = i + random.nextInt(length-i);
            int temp = array[i];
            array[i] = array[ran];
            array[ran] = temp;
        }
        return array;
    }

    // Shingle calculation
    @Override
    public int[] call() throws Exception {
        int[] ans = new int[numNodes];
        int[] superAns = new int[numSuperNodes];
        Int2IntOpenHashMap hash = new Int2IntOpenHashMap();

        // Do random Permute Here
        int[] perm = getRandomPermutation(numNodes);
        int tmp;
        // loop over all the nodes, calculate shingles for them
        for(int i=0;i<numNodes;i++) ans[i] = perm[i];
        for(long e: edges) {
            if(e < 0) continue;
            int src = (int)(e >> 32);
            int dst = (int)(e & 0x7FFFFFFFL);
            ans[src] = (perm[dst] < ans[src]) ? perm[dst] : ans[src];
            ans[dst] = (perm[src] < ans[dst]) ? perm[src] : ans[dst];
        }

        for(int i=0;i<numSuperNodes;i++) superAns[i] = 0x7FFFFFFF;
        for(int i=0;i<numSuperNodes;i++){
            int supernode = superList.getInt(i);
            superAns[i] = Integer.MAX_VALUE;
            for(int j: insideSupernode[supernode]){
                superAns[i] = (ans[j] < superAns[i]) ? ans[j] : superAns[i];
            }
        }
        return superAns;
    }
}
