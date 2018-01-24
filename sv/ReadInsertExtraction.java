import java.util.*;
import java.io.*;
public class ReadInsertExtraction {
public static void main(String[] args) throws IOException
{
	String readFn = args[0];
	int pos = Integer.parseInt(args[1]);
	Scanner input = new Scanner(new FileInputStream(new File(readFn)));
	boolean whole = false;
	if(args.length > 2 && args[2].equals("whole")) whole = true;
	PrintWriter out = new PrintWriter(System.out);
	
	int minLength = 30;
	int maxDist = 500;
	int flank = 1000;
	
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.charAt(0) == '@') continue;
		StringTokenizer str = new StringTokenizer(line);
		ArrayList<String> tokens = new ArrayList<String>();
		while(str.hasMoreTokens())
		{
			tokens.add(str.nextToken());
		}
		String name = tokens.get(0);
		int refIdx = Integer.parseInt(tokens.get(3));
		String cigar = tokens.get(5);
		String seq = tokens.get(9);
		parse(cigar);
		int idx = 0;
		int n = alignmentTypes.length;
		int startPos = -1, endPos = -1;
		int bestDist = maxDist+1;
		boolean softClipped = false;
		for(int i = 0; i<n; i++)
		{
			int len = alignmentLengths[i];
			char t = alignmentTypes[i];
			if(t == 'I' || t == 'S')
			{
				//if(len > 10) System.out.println(len+" "+refIdx);
				if(len >= minLength && Math.abs(refIdx - pos) <= Math.min(maxDist, bestDist))
				{
					bestDist = Math.abs(refIdx - pos);
					startPos = idx;
					endPos = idx + len;
					if(t == 'S') softClipped = true;
					else softClipped = false;
				}
			}
			if(t != 'D' && t != 'H') idx += len;
			if(t != 'I' && t != 'S' && t != 'H') refIdx += len;
		}
		if(startPos != -1)
		{
			int startIdx = Math.max(0, startPos - flank);
			int endIdx = Math.min(idx, endPos + flank);
			System.out.println(">" + name);
			System.out.println(seq.substring(startIdx, endIdx));
			//System.out.println((whole || softClipped) ? seq : seq.substring(startIdx, endIdx));
		}
	}
	
	out.close();
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
