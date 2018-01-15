import java.util.*;
import java.io.*;
public class CoordsParser {
public static void main(String[] args) throws IOException
{
	Scanner stdIn = new Scanner(System.in);
	int pos = Integer.parseInt(args[2]);
        String fn = args[0];
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	ArrayList<Match> list = new ArrayList<Match>();
	String fastaFile = args[1];
	Scanner fastaInput = new Scanner(new FileInputStream(new File(fastaFile)));
	fastaInput.nextLine();
	String seq = fastaInput.next();
	PrintWriter logFile = new PrintWriter(new File(fastaFile+".log"));
	while(input.hasNext())
	{
		String s = input.nextLine();
		if(s.startsWith("=")) break;
	}
	while(input.hasNextLine())
	{
		int s1 = input.nextInt();
		int e1 = input.nextInt();
		input.next();
		int s2 = input.nextInt();
		int e2 = input.nextInt();
		input.nextLine();
		list.add(new Match(s1, e1, s2, e2));
		logFile.println("new match: " + s1+" "+e1+" "+s2+" "+e2);
	}
	int best = -1, best2 = -1;
	int other = -1, other2 = -1;
	int p = -1, p2 = -1;
	for(int i = 0; i<list.size(); i++)
	{
		Match cur = list.get(i);
		int sdist = Math.abs(cur.s1 - pos);
		int edist = Math.abs(cur.e1 - pos);
		boolean startBetter = sdist < edist;
		if(startBetter)
		{
			if(best == -1 || sdist < best)
			{
				best2 = best;
				other2 = other;
				p2 = p;
				best = sdist;
				other = cur.s2;
				p = cur.e1;
			}
			else if(best2 == -1 || sdist < best2)
			{
				best2 = sdist;
				other2 = cur.s2;
				p2 = cur.e1;
			}
		}
		else
		{
			if(best == -1 || edist < best)
			{
				best2 = best;
				other2 = other;
				p2 = p;
				best = edist;
				other = cur.e2;
				p = cur.s1;
			}
			else if(best2 == -1 || edist < best2)
			{
				best2 = edist;
				other2 = cur.e2;
				p2 = cur.s1;
			}
		}
	}
	logFile.println("other: " + other + ", other2: " + other2);
	logFile.println("p: " + p + ", p2: " + p2); 
	if(other == -1 || other2 == -1) return;
	String s = seq.substring(Math.min(other2, other), Math.max(other2, other));
	if(p > p2)
	{
		int temp = other;
		other = other2;
		other2 = temp;
	}
	// Now the substring from other to other2 is what we want
	String res = "";
	System.out.println(res = (other < other2 ? s : reverse(s)));
	logFile.println(seq.contains(res));
	logFile.close();
}
static String reverse(String s)
{
	s = s.replaceAll("A", "#");
	s = s.replaceAll("T", "A");
	s = s.replaceAll("#", "T");
	s = s.replaceAll("C", "#");
	s = s.replaceAll("G", "C");
	s = s.replaceAll("#", "G");
	return new StringBuilder(s).reverse().toString();
}
static class Match
{
	int s1, e1;
	int s2, e2;
	Match(int ss1, int ee1, int ss2, int ee2)
	{
		s1 = ss1; e1 = ee1; s2 = ss2; e2 = ee2;
	}
}
}
