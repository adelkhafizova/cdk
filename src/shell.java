import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;









import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemModel;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.signature.AtomSignature;

class shell {
	public void testSDFFile() throws FileNotFoundException {
        String filename = "C:\\Users\\Vitaliy\\Downloads\\Ligands_noHydrogens_withMissing_1_Instances.sdf"; // a multi molecule SDF file
        InputStream ins = new FileInputStream(filename);
        try {
        	 File sdfFile = new File("C:\\Users\\Vitaliy\\Downloads\\0008,1,0,0013.sdf");
        	 IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(sdfFile), DefaultChemObjectBuilder.getInstance());
        	 while (reader.hasNext()) {
        		 IMolecule molecule = (IMolecule)reader.next();
        		 System.out.println("More molecules");
        		 for (int i = 0; i < molecule.getAtomCount(); i++) {
        			 AtomSignature aa = new AtomSignature(i, 1, molecule);
           			 String aa_s = new String(aa.toString());
        			 System.out.println(aa_s);
        		 }
        	}
            /*MDLReader reader = new MDLReader(ins);
            ChemFile fileContents = (ChemFile)reader.read(new ChemFile());
            if (1 == fileContents.getChemSequenceCount()) {
            	System.out.println("OK");
            }
            else
            	System.out.println("Not OK");
            org.openscience.cdk.interfaces.IChemSequence sequence = fileContents.getChemSequence(0);
            if (sequence == null) {
            	System.out.println("Not OK");
            }*/
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
    }

	public static void main(String[] args) throws FileNotFoundException {
		shell e = new shell();
		e.testSDFFile();
	}

}
