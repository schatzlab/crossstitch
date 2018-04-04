import java.util.*;
import java.io.*;
public class VCFEditor {
public static void main(String[] args) throws IOException
{
	String readFn = args[0], vcfFn = args[1];
	String outputFn = args[2];
	Scanner readsInput = new Scanner(new FileInputStream(new File(readFn)));
	Scanner vcfScanner = new Scanner(new FileInputStream(new File(vcfFn)));
	PrintWriter out = new PrintWriter(new File(outputFn));
	PrintWriter logout = new PrintWriter(new File(vcfFn+".log"));
	HashMap<String, String> readMap = new HashMap<String, String>();
	String find = "<INS>";
	while(readsInput.hasNext())
	{
		String nextLine = readsInput.nextLine();
		while(readsInput.hasNext() && nextLine.charAt(0) != '>')
			nextLine = readsInput.nextLine();
		String id = nextLine.substring(1);
		
		String seq = readsInput.nextLine();
		//System.out.println(id+" "+seq);
		//if(id.charAt(0) != 'i') return;
		readMap.put(id.substring(id.indexOf('_')+1), seq);
	}
	while(vcfScanner.hasNext())
	{
		String line = vcfScanner.nextLine();
		if(line.charAt(0) == '#')
		{
		    out.println(line);
		    continue;
		}
		String[] tokens = line.split("\t");
		int index = 0;
		int curPos = -1;
		String ch = "";
		for(String t : tokens)
		{
		    if(t.contains("chr"))
		    {
			ch = t;
		        curPos = Integer.parseInt(tokens[index+1]);
		        break;
		    }
			index++;
		}
		if(curPos == -1 && tokens.length > 1)
		{
			ch = tokens[0];
			curPos = Integer.parseInt(tokens[1]);
		}
		if(curPos == -1) continue;
		if(!line.contains(find))
		{
		    out.println(line);
		    continue;
		}
		if(!readMap.containsKey(ch+":"+curPos+""))
		{
			logout.println(ch+":"+curPos+"");
			continue;
		}
		System.out.println(curPos);
		String ins = readMap.get(ch+":"+curPos+"");
		if(ins.equals("null"))
		{
			logout.println(ch+":"+curPos+" "+ins);
			continue;
		}
		for(int i = 0; i<line.length()-4; i++)
		{
			if(line.substring(i, i+4).equals("SEQ="))
			{
				String before = line.substring(0,i+4);
				String after = line.substring(i+4);
				if(ins != null && ins.length() > 0)
				{
				    after = readMap.get(ch+":"+curPos+"") + after.substring(after.indexOf(';'));
				}
				else
				{
				    logout.println("Used sniffles output for " + (ch+":"+curPos+""));
				}
				out.println(before+after);
				break;
			}
		}
	}
	
	out.close();
	logout.close();
}
}
