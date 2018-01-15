import java.util.*;
import java.io.*;
public class BestInsertFinder {
public static void main(String[] args) throws IOException
{
	int pos = Integer.parseInt(args[0]);
	int n = args.length - 1;
	if(n == 0) return;
	String[] fns = new String[n];
	for(int i = 0; i<n; i++) fns[i] = args[i+1];
	String[] seqs = new String[n];
	int besti = 0;
	for(int i = 0; i<n; i++)
	{
		if(!fns[i].contains(pos+"_")) continue;
		String cur = fns[i];
		Scanner input = new Scanner(new FileInputStream(new File(cur)));
		input.nextLine();
		seqs[i] = input.nextLine();
		if(seqs[besti] == null || seqs[i].length() > seqs[besti].length())
		{
			besti = i;
		}
	}
	System.out.println(seqs[besti]);
}
}
