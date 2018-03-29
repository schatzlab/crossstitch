import java.util.*;
import java.io.*;
public class BestInsertFinder2 {
public static void main(String[] args) throws IOException
{
	String samFn = args[0];
	int pos = Integer.parseInt(args[1]);
	int offset = Integer.parseInt(args[2]);
	String type = args[3]; // SEQ or POS
	int minLength = args.length > 4 ? Integer.parseInt(args[4]) : 30;
	int maxDist = args.length > 5 ? Integer.parseInt(args[5]) : 5000;
	Scanner input = new Scanner(new FileInputStream(new File(samFn)));
	String insertSeq = "";
	String bestTig = "";
	int middle = 0;
	int bestScore = -987654321;
	int genomeLength = -1;
	int bestPos = pos;
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.charAt(0) == '@')
		{
			if(line.substring(0, 3).equals("@SQ"))
			{
				StringTokenizer str = new StringTokenizer(line);
				while(str.hasMoreTokens())
				{
					String cur = str.nextToken();
					if(cur.length() >= 3 && cur.substring(0, 3).equals("LN:"))
					{
						genomeLength = Integer.parseInt(cur.substring(3));
						middle = pos < offset ? pos+1 : (offset + 1);
						
					}
				}
			}
			continue;
		}
		StringTokenizer str = new StringTokenizer(line);
		ArrayList<String> tokens = new ArrayList<String>();
		while(str.hasMoreTokens())
		{
			tokens.add(str.nextToken());
		}
		String tigName = tokens.get(0);
		String cigar = tokens.get(5);
		String seq = tokens.get(9);
		int refStartIdx = Integer.parseInt(tokens.get(3))-1;
		parse(cigar);
		int idx = 0;
		int refIdx = refStartIdx;
		int n = alignmentTypes.length;
		for(int i = 0; i<n; i++)
		{
			int len = alignmentLengths[i];
			char t = alignmentTypes[i];
			if(t == 'I')
			{
				if(score(refIdx, middle, len) > bestScore && len >= minLength && Math.abs(refIdx - middle) <= maxDist)
				{
				    bestPos = refIdx;
					bestScore = score(refIdx, middle, len);
					insertSeq = seq.substring(idx, idx+len);
					bestTig = tigName;
				}
			}
			if(t != 'D' && t != 'H' && t != 'N' && t != 'P') idx += len;
			if(t != 'I' && t != 'S' && t != 'H' && t != 'P') refIdx += len;
		}
	}
	if(type.equals("SEQ")) System.out.println(insertSeq);
	else System.out.println(bestPos);
}
static int score(int pos, int middle, int len)
{
	return 20*len - Math.abs(pos-middle)*Math.abs(pos-middle);
}
static char[] alignmentTypes;
static int[] alignmentLengths;
static void parse(String cigar)
{
	ArrayList<Integer> lens = new ArrayList<Integer>();
	ArrayList<Character> types = new ArrayList<Character>();
	
	int idx = 0;
	if(cigar.length() > 1)
	{
		while(idx < cigar.length())
		{
			int nidx = idx;
			while(cigar.charAt(nidx) <= '9' && cigar.charAt(nidx) >= '0') nidx++;
			lens.add(Integer.parseInt(cigar.substring(idx, nidx)));
			types.add(cigar.charAt(nidx));
			idx = nidx+1;
		}
	}
	
	alignmentTypes = new char[types.size()];
	for(int i = 0; i<types.size(); i++) alignmentTypes[i] = types.get(i);
	
	alignmentLengths = new int[lens.size()];
	for(int i = 0; i<lens.size(); i++) alignmentLengths[i] = lens.get(i);
}
}
