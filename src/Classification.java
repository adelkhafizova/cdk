import java.io.*;
import java.util.*;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.signature.AtomSignature;
import org.openscience.cdk.signature.MoleculeSignature;
import org.openscience.cdk.smiles.SmilesGenerator;

class Classification {
    Classification() {
        this.molecule_data = new ArrayList<IMolecule>();
        this.heights_active = new HashMap<String, Integer>();
        this.heights_inactive = new HashMap<String, Integer>();
        signatures_active = new HashMap<String, String>();
        signatures_inactive = new HashMap<String, String>();
    }

    void read_signatures(String path, HashMap<String, String> signatures) throws IOException {
        InputStream fis = new FileInputStream(path);
        String pair;
        InputStreamReader isr = new InputStreamReader(fis);
        BufferedReader br = new BufferedReader(isr);
        while ((pair = br.readLine()) != null) {
            String[] parts = pair.split(" ");
            signatures.put(parts[0], parts[1]);
        }
        br.close();
        isr.close();
        fis.close();
    }

    void read_signatures_heights(String path, HashMap<String, Integer> signatures) throws IOException {
        InputStream fis = new FileInputStream(path);
        String pair;
        InputStreamReader isr = new InputStreamReader(fis);
        BufferedReader br = new BufferedReader(isr);
        while ((pair = br.readLine()) != null) {
            String[] parts = pair.split(" ");
            signatures.put(parts[0], Integer.parseInt(parts[1]));
        }
        br.close();
        isr.close();
        fis.close();
    }

    void read_molecules(String path, ArrayList<IMolecule> molecule_data) throws IOException {
        File sdfFile = new File(path);
        IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(sdfFile),
                                                           DefaultChemObjectBuilder.getInstance());
        while (reader.hasNext()) {
            IMolecule molecule = (IMolecule) reader.next();
            molecule_data.add(molecule);
        }
    }

    void data_initialization(String path) throws IOException {
        read_signatures(path + "final_signatures_active.txt", signatures_active);
        read_signatures(path + "final_signatures_inactive.txt", signatures_inactive);
        read_signatures_heights(path + "signatures_active_heights.txt", heights_active);
        read_signatures_heights(path + "signatures_inactive_heights.txt", heights_inactive);
    }

    String atom_signature_to_inchi(AtomSignature as) throws CDKException {
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        IAtomContainer ac = MoleculeSignature.fromSignatureString(as.toCanonicalString(), builder);
        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
        InChIGenerator gen = factory.getInChIGenerator(ac);
        String inchi = gen.getInchi();
        return inchi;
    }

	void classify_data(String molecule_file) throws CDKException, IOException {
		ArrayList<Boolean> active = new ArrayList<Boolean>();
		ArrayList<Boolean> inactive = new ArrayList<Boolean>();
		ArrayList<HashSet<String>> active_fragments = new ArrayList<HashSet<String>>();
		ArrayList<HashSet<String> > inactive_fragments = new ArrayList<HashSet<String>>();
        read_molecules(new String(molecule_file), molecule_data);

		int index = 0;
		System.out.println("Progress output");
		for (IMolecule m : molecule_data) {
			active.add(false);
			inactive.add(false);
			active_fragments.add(new HashSet<String>());
			inactive_fragments.add(new HashSet<String>());
			for (int i = 0; i < m.getAtomCount(); ++i) {
				for (Map.Entry<String, String> substructure: signatures_active.entrySet()) {
					String as = substructure.getKey();
					int begin = as.indexOf('[');
					int end = as.indexOf(']', begin + 1);
					if (m.getAtom(i).getSymbol().equals(as.substring(begin + 1, end))) {
                        AtomSignature atom_s = new AtomSignature(i, heights_active.get(as), m);
						String InChi = atom_signature_to_inchi(atom_s);
						if (InChi.equals(substructure.getValue())) {
							active_fragments.get(index).add(InChi);
							active.set(index, true);
						}
					}
				}
				for (Map.Entry<String, String> substructure: signatures_inactive.entrySet()) {
					String as = substructure.getKey();
					int begin = as.indexOf('[');
					int end = as.indexOf(']', begin + 1);
					if (m.getAtom(i).getSymbol().equals(as.substring(begin + 1, end))) {
                        AtomSignature atom_s = new AtomSignature(i, heights_inactive.get(as), m);
                        String InChi = atom_signature_to_inchi(atom_s);
						if (InChi.equals(substructure.getValue())) {
							inactive_fragments.get(index).add(InChi);
							inactive.set(index, true);
						}
					}
				}
			}
			++index;
			System.out.println(((double) index) / this.molecule_data.size() * 100);
		}
		int active_count = 0;
		int inactive_count = 0;
		ArrayList<Integer> final_activity = new ArrayList<Integer>();
		for (int i = 0; i < active.size(); ++i) {
			if (active.get(i)) {
				if (inactive.get(i)) {
					System.out.println("unrecognized");
					final_activity.add(2);
				} else {
					System.out.println("mutagen");
					++active_count;
					final_activity.add(1);
				}
			} else {
				if (inactive.get(i)) {
					System.out.println("nonmutagen");
					++inactive_count;
					final_activity.add(0);
				} else {
					System.out.println("unrecognized");
					final_activity.add(2);
				}
			}
		}
		System.out.println("Active count" + active_count);
		System.out.println("Inactive count" + inactive_count);
		PrintWriter output_file = new PrintWriter("final_results_new", "UTF-8");
		SmilesGenerator generator=new SmilesGenerator();
		for (int i = 0; i < molecule_data.size(); ++i) {
			output_file.println(generator.createSMILES(molecule_data.get(i)));
			//System.out.println(generator.createSMILES(molecule_data.get(i)));
			if (final_activity.get(i) == 0) {
				output_file.println("nonmutagen" + '\n');
				//System.out.println("nonmutagen" + '\n');
			} else {
				if (final_activity.get(i) == 1)  {
					output_file.println("mutagen" + '\n');
					//System.out.println("mutagen" + '\n');
				} else {
					output_file.println("unrecognized" + '\n');
					//System.out.println("unrecognized" + '\n');
				}
			}
			Iterator<String> itr = active_fragments.get(i).iterator();
			output_file.println("Active fragments");
			while (itr.hasNext()) {
				output_file.println(itr.next() + ' ');
				//System.out.println(itr.next() + ' ');
			}
			output_file.println();
			output_file.println("Inactive fragments");
			itr = inactive_fragments.get(i).iterator();
			while (itr.hasNext()) {
				output_file.println(itr.next() + ' ');
				//System.out.println(itr.next() + ' ');
			}
			output_file.println();
		}
		output_file.close();
	}
	ArrayList<IMolecule> molecule_data;
    HashMap<String, Integer> heights_active;
    HashMap<String, Integer> heights_inactive;
    HashMap<String, String> signatures_active;
    HashMap<String, String> signatures_inactive;
}
