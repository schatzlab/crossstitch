import java.util.*;
import java.io.*;
public class FastaFileFixer2 {
public static void main(String[] args) throws IOException
{
	Scanner input = new Scanner(new FileInputStream(new File(args[0])));
	PrintWriter out = new PrintWriter(new File(args[1]));
	String name = "";
	StringBuilder seq = new StringBuilder("");
	while(input.hasNext())
	{
		String s = input.nextLine();
		if(s.startsWith(">"))
		{
			if(name.length()>0)
			{
				out.println(">" + name);
				out.println(seq);
				seq = new StringBuilder();
			}
			name = s.substring(1);
		}
		else seq.append(s);
	}
	if(name.length()>0)
	{
	    out.println(">" + name);
		out.println(seq);
	}
		
		
	
	
	out.close();
}
}
