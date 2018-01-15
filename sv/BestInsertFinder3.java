import java.util.*;
import java.io.*;
public class BestInsertFinder3 {
public static void main(String[] args) throws IOException
{
	String fn = args.length > 0 ? args[0] : "/home/mkirsche/study/sniffles_100kbp.vcf";
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	int bestDist = (int)2e9;
	String insert = "";
	int maxDist = 500;
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.charAt(0) == '#') continue;
		if(!line.contains("<INS>") || !line.contains("SEQ=")) continue;
		StringTokenizer str = new StringTokenizer(line);
		ArrayList<String> tokens = new ArrayList<String>();
		while(str.hasMoreTokens())
		{
			tokens.add(str.nextToken());
		}
		String chr = tokens.get(0);
		int mid = getMid(chr);
		int pos = Integer.parseInt(tokens.get(1));
		//System.out.println(mid+" "+pos);
		int dist = Math.abs(mid - pos);
		if(dist < bestDist && dist <= maxDist)
		{
			bestDist = dist;
			String insertSuffix = line.substring(line.indexOf("SEQ=")+4);
			insert = insertSuffix.substring(0, insertSuffix.indexOf(';'));
		}
	}
	System.out.println(insert);
}
static int getMid(String chr)
{
	String range = chr.substring(chr.indexOf(':')+1);
	String start = range.substring(0,range.indexOf('-'));
	String end = range.substring(range.indexOf('-')+1);
	int s = Integer.parseInt(start);
	int e = Integer.parseInt(end);
	int mid = (e - s)/2;
	return mid;
}
}
