import java.util.*;
import java.io.*;
public class VCFEditor {
public static void main(String[] args) throws IOException
{
	String readFn = args[0], posFn = args[1], vcfFn = args[2];
	String outputFn = args[4];
	String fastaFn = args[3];
	int left = Integer.parseInt(args[5]);
	int right = Integer.parseInt(args[6]);
	Scanner readsInput = new Scanner(new FileInputStream(new File(readFn)));
	Scanner vcfScanner = new Scanner(new FileInputStream(new File(vcfFn)));
	Scanner posInput = new Scanner(new FileInputStream(new File(posFn)));
	PrintWriter out = new PrintWriter(new File(outputFn));
	PrintWriter logout = new PrintWriter(new File(vcfFn+".log"));
	HashMap<String, String> readMap = new HashMap<String, String>();
	String find = "<INS>";
	while(readsInput.hasNext())
	{
		String id = readsInput.nextLine().substring(1);
		String seq = readsInput.nextLine();
		readMap.put(id.substring(id.indexOf('_')+1), seq);
	}
	HashMap<String, Integer> posMap = new HashMap<String, Integer>();
	while(posInput.hasNext())
	{
		String id = posInput.nextLine().substring(1);
		int p = Integer.parseInt(posInput.nextLine());
		posMap.put(id.substring(id.indexOf('_')+1), p);
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
		int curPos = Integer.parseInt(tokens[index+1]);
		String ch = tokens[0];
		if(!line.contains(find))
		{
		    if(line.contains("<DEL>"))
		    {
		        String seq = "";
		        if(line.contains("SEQ=")) seq = getField(line, "SEQ");
		        else
		        {
		            int start = curPos;
		            int end = start + Integer.parseInt(getField(line, "SVLEN"));
		            String command = "samtools faidx " + fastaFn + " " + ch + ":" + (start) + "-" + (end-1);
		            Process child = Runtime.getRuntime().exec(command);
                    InputStream seqStream = child.getInputStream();
                    Scanner seqInput = new Scanner(seqStream);
                    seqInput.next();
                    while(seqInput.hasNext()) seq += seqInput.next();
		        }
		        StringTokenizer str = new StringTokenizer(line);
		        int tokenIdx = 0;
		        while(str.hasMoreTokens())
		        {
		            String cur = str.nextToken();
	                if(tokenIdx == 3)
	                {
	                    cur = seq;
	                }
	                else if(tokenIdx == 4)
	                {
	                    cur = "X";
	                }
		            out.print(cur);
		            if(str.hasMoreTokens()) out.print('\t');
		            else out.println();
		            tokenIdx++;
		        }
		    }
		    else out.println(line);
		    continue;
		}
		if(!readMap.containsKey(ch+":"+curPos+""))
		{
			logout.println(ch+":"+curPos+"");
			continue;
		}
		String ins = readMap.get(ch+":"+curPos+"");
		Integer svPos = posMap.get(ch+":"+curPos+"");
		if(ins == "null")
		{
			logout.println(ch+":"+curPos+" "+ins);
			continue;
		}
		String output = line;
		if(ins.length() == 0)
		{
		    logout.println("Used sniffles output for " + (ch+":"+curPos+""));
		}
		else
		{
		    output = substitute(line, "SEQ", ins);
		    output = substitute(output, "SVLEN", ins.length()+"");
		}
		StringTokenizer str = new StringTokenizer(output);
		int tokenIdx = 0;
		String command = "samtools faidx " + fastaFn + " " + ch + ":" + (svPos-left+1) + "-" + (svPos+right);
	System.out.println(command);
	Process child = Runtime.getRuntime().exec(command);
        InputStream seqStream = child.getInputStream();
        Scanner seqInput = new Scanner(seqStream);
        seqInput.next();
        String seq = seqInput.next();
		while(str.hasMoreTokens())
		{
		    String cur = str.nextToken();
		    if(ins.length() > 0)
		    {
		        if(svPos != null && svPos != -1 && tokenIdx == 1)
		        {
		            cur = svPos+"";
		        }
		        if(tokenIdx == 3)
		        {
		            cur = seq;
		        }
		        else if(tokenIdx == 4)
		        {
		            cur = seq.substring(0, left) + ins + seq.substring(left);
		        }
		    }
		    out.print(cur);
		    if(str.hasMoreTokens()) out.print('\t');
		    else out.println();
		    tokenIdx++;
		}
	}
	
	out.close();
	logout.close();
}
static String getField(String s, String field)
{
    int n = s.length(), m = field.length();
    for(int i = 0; i+m+1 <= n; i++)
    {
        if(s.substring(i, i+m+1).equals(field+"="))
        {
            String after = s.substring(i+m+1);
            return after.substring(0, after.indexOf(';'));
        }
    }
    return "";
}
static String substitute(String s, String field, String val)
{
    if(val.length() == 0)
    {
        return s;
    }
    int n = s.length(), m = field.length();
    for(int i = 0; i+m+1 <= n; i++)
    {
        if(s.substring(i, i+m+1).equals(field+"="))
        {
            String before = s.substring(0, i+m+1);
            String after = s.substring(i+m+1);
            after = after.substring(after.indexOf(';'));
            return before + val + after;
        }
    }
    return "";
}
}
