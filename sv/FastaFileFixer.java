import java.util.*;
import java.io.*;
public class FastaFileFixer {
public static void main(String[] args) throws IOException
{
	Scanner input = new Scanner(new FileInputStream(new File(args[0])));
	PrintWriter out = new PrintWriter(new File(args[1]));
	while(input.hasNext())
	{
		String s = input.nextLine().substring(1);
		String t = input.nextLine();
		input.nextLine();
		input.nextLine();
		out.println(">"+s+"\n"+t);
	}
	out.close();
}
}
