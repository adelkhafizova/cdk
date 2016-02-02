import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.io.*;
//import java.math.BigInteger;
import java.util.*;
import java.util.List;
import java.util.concurrent.ExecutionException;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.exception.IncorrectUseOfCDKCoreClassError;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.formats.SybylDescriptorFormat;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.renderer.AtomContainerRenderer;
import org.openscience.cdk.renderer.RendererModel;
import org.openscience.cdk.renderer.font.AWTFontManager;
import org.openscience.cdk.renderer.generators.BasicAtomGenerator;
import org.openscience.cdk.renderer.generators.BasicBondGenerator;
import org.openscience.cdk.renderer.generators.BasicSceneGenerator;
import org.openscience.cdk.renderer.generators.IGenerator;
import org.openscience.cdk.renderer.visitor.AWTDrawVisitor;
import org.openscience.cdk.signature.AtomSignature;
import org.openscience.cdk.signature.MoleculeSignature;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.io.SDFWriter;

import javax.imageio.ImageIO;
//import com.google.common.math.BigIntegerMath;

class p_value_pair {
    public p_value_pair(Double p_value_sign, Double p_value_int) {
        this.p_value_first = p_value_sign;
        this.p_value_second = p_value_int;
    }
    Double p_value_first;
    Double p_value_second;
    public String toString() {
        String local = new String(p_value_first.toString() + ' ' + p_value_second.toString());
        return local;
    }
}


public class Algorithm {
    Algorithm(double p_value_threshold, int minimum_occurrence, double substructure_frequency,
              HashMap<String, Double> ex_active, HashMap<String, Double> ex_inactive) {
        this.p_value_threshold = p_value_threshold;
        this.minimum_occurrence = minimum_occurrence;
        this.substructure_frequency = substructure_frequency;
        this.molecule_data = new ArrayList<IMolecule>();
        this.activity_data = new ArrayList<Boolean>();
        this.atom_states = new ArrayList<ArrayList<Boolean>>(); //true - nothing else to do with it, false - look
        this.pi = 3.1415926535897932384;
        this.final_signatures_active = new HashMap<String, Double>();
        this.final_signatures_inactive = new HashMap<String, Double>();
        this.active = new HashMap<String, String>();
        this.inactive = new HashMap<String, String>();
        this.heights_active = new HashMap<String, Integer>();
        this.heights_inactive = new HashMap<String, Integer>();
        this.data_num = 0;
        this.active_num = 0;
        this.prev_active = ex_active;
        this.prev_inactive = ex_inactive;
        this.int_active = new HashMap<String, p_value_pair >();
        this.int_inactive = new HashMap<String, p_value_pair>();
    }

    void data_initialization(String pathname, String activity_name, String inactivity_name) {
        try {
            File sdfFile = new File(pathname);
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(sdfFile), DefaultChemObjectBuilder.getInstance());
            int inactive_num = 0;
            while (reader.hasNext()) {
                IMolecule molecule = (IMolecule) reader.next();
                String properties = molecule.getProperties().values().toString();
                if (properties.contains(inactivity_name)) {
                    activity_data.add(false);
                    ++inactive_num;
                    molecule_data.add(molecule);
                    ArrayList<Boolean> atom_states_local = new ArrayList<Boolean>(molecule.getAtomCount());
                    for (int t = 0; t < molecule.getAtomCount(); t++) {
                        atom_states_local.add(false);
                    }
                    atom_states.add(atom_states_local);
                    data_num++;
                } else {
                    if (properties.contains(activity_name)) {
                        activity_data.add(true);
                        ++active_num;
                        molecule_data.add(molecule);
                        ArrayList<Boolean> atom_states_local = new ArrayList<Boolean>(molecule.getAtomCount());
                        for (int t = 0; t < molecule.getAtomCount(); t++) {
                            atom_states_local.add(false);
                        }
                        atom_states.add(atom_states_local);
                        data_num++;
                    }
                }
            }
            if (active_num + inactive_num != molecule_data.size()) {
                System.err.println(active_num);
                System.err.println(inactive_num);
                System.err.println(molecule_data.size());
                throw new Exception("Invalid data, not all molecules have activity value");
            }
            //here fill molecule_data, activity_data and allocate memory for atom states and fill it with 0
            //while initializing determine how many activity compounds there are
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
    }

    class Signature_state {
        Signature_state() {
            this.met_positions = new ArrayList<Signature_position>();
            n_dash = 0;
            m_dash = 0;
        }

        void add(int i, int j, Boolean activity) {
            Signature_position sp = new Signature_position(i, j);
            this.met_positions.add(sp);
            n_dash++;
            if (activity)
                m_dash++;
        }

        //int processed_condition; //not met; met but not interesting; met but interesting; included in final list
        ArrayList<Signature_position> met_positions; //positions of atoms, where this signature was met
        double n_dash;//number of substructures
        double m_dash;//number of active substructures
    }

    class Signature_position {
        int molecule_number;
        int atom_number;

        Signature_position(int i, int j) {
            this.molecule_number = i;
            this.atom_number = j;
        }
    }

    String atom_signature_to_inchi(AtomSignature as) throws CDKException {
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        IAtomContainer ac = MoleculeSignature.fromSignatureString(as.toCanonicalString(), builder);
        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
        InChIGenerator gen = factory.getInChIGenerator(ac);
        return gen.getInchi();
    }

    String atom_signature_to_inchi(String as) throws CDKException {
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        IAtomContainer ac = MoleculeSignature.fromSignatureString(as, builder);
        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
        InChIGenerator gen = factory.getInChIGenerator(ac);
        return gen.getInchi();
    }

    void find_substructures_with_height(Map<String, Signature_state> current_signatures,
                                        Map<String, String> signature_to_inchi, int height) throws CDKException {
        int molecule_index = 0;
        for (IMolecule m : this.molecule_data) { //find all interesting substructures
            for (int i = 0; i < m.getAtomCount(); i++) {
                if (this.atom_states.get(molecule_index).get(i) == false) { //not in the current signatures and not viewed
                    AtomSignature as = new AtomSignature(i, height, m);
                    signature_to_inchi.put(as.toCanonicalString(), atom_signature_to_inchi(as));
                    if (current_signatures.containsKey(as.toString()) == false) {
                        Signature_state ss = new Signature_state();
                        ss.add(molecule_index, i, this.activity_data.get(molecule_index));
                        current_signatures.put(as.toString(), ss);
                    } else {
                        current_signatures.get(as.toString()).add(molecule_index, i, this.activity_data.get(molecule_index));
                    }
                }
            }
            molecule_index++;
        }
    }

    void process_signatures(Map<String, Signature_state> current_signatures, int height) throws CDKException, FileNotFoundException, UnsupportedEncodingException {
        for (Map.Entry<String, Signature_state> substructure : current_signatures.entrySet()) {
            if (substructure.getValue().n_dash > minimum_occurrence) {
                double current_frequency = substructure.getValue().m_dash / substructure.getValue().n_dash;
                double n_ = substructure.getValue().n_dash;
                double m_ = substructure.getValue().m_dash;
                if (current_frequency > substructure_frequency) {
                    double current_p_value_active = 0;
                    for (double k = m_; k < n_; k++) {
                        if (n_ - k == 0 || k == 0 || data_num == 0) {
                            System.out.println(n_ - k + ' ' + k + ' ' + data_num);
                        }
                        current_p_value_active += Math.pow(n_ / (2 * pi * (n_ - k) * k), 0.5) * Math.pow(active_num / k, k) *
                                Math.pow((data_num - active_num) / (n_ - k), n_ - k) * Math.pow(n_ / data_num, n_);
                        //BigInteger coef = BigIntegerMath.binomial((int) n_, (int) k);
                        //true_p_value += (coef.doubleValue() * Math.pow(active_num, k)) / Math.pow(data_num, k) *
                        //        Math.pow((data_num - active_num), n_ - k) / Math.pow(data_num, n_ - k);
                    }
                    current_p_value_active += Math.pow(active_num / data_num, n_);
                    if (current_p_value_active < p_value_threshold) {
                        for (Signature_position sp : substructure.getValue().met_positions) {
                            atom_states.get(sp.molecule_number).set(sp.atom_number, true);
                        }

                        String as_string = substructure.getKey();
                        String InChi = atom_signature_to_inchi(as_string);
                        final_signatures_active.put(InChi, current_p_value_active);
                        active.put(InChi, as_string);
                        heights_active.put(as_string, height);
                    }
                    else {
                        if (prev_active.size() != 0) {
                            String as_string = substructure.getKey();
                            String InChi = atom_signature_to_inchi(as_string);
                            if (prev_active.containsKey(InChi)) {
                                int_active.put(InChi, new p_value_pair(prev_active.get(InChi), current_p_value_active));
                        }
                        }
                    }
                }

                current_frequency = (substructure.getValue().n_dash - substructure.getValue().m_dash) / substructure.getValue().n_dash;

                if (current_frequency > substructure_frequency) {
                    double current_p_value_inactive = 0;
                    m_ = n_ - substructure.getValue().m_dash;
                    for (double k = m_; k < n_; k++) {
                        current_p_value_inactive += Math.pow(n_ / (2 * pi * (n_ - k) * k), 0.5) * Math.pow((data_num - active_num) / k, k) *
                                Math.pow(active_num / (n_ - k), n_ - k) * Math.pow(n_ / data_num, n_);
                    }
                    current_p_value_inactive += Math.pow((data_num - active_num) / data_num, n_);
                    if (current_p_value_inactive < p_value_threshold) {
                        for (Signature_position sp : substructure.getValue().met_positions) {
                            atom_states.get(sp.molecule_number).set(sp.atom_number, true);
                        }
                        String as_string = substructure.getKey();
                        String InChi = atom_signature_to_inchi(as_string);
                        final_signatures_inactive.put(InChi, current_p_value_inactive);
                        inactive.put(InChi, as_string);
                        heights_inactive.put(as_string, height);
                    }
                    else {
                        if (prev_inactive.size() != 0) {
                            String as_string = substructure.getKey();
                            String InChi = atom_signature_to_inchi(as_string);
                            if (prev_inactive.containsKey(InChi)) {
                                int_inactive.put(InChi, new p_value_pair(prev_inactive.get(InChi), current_p_value_inactive));
                            }
                        }
                    }
                }
            } else {
                for (Signature_position sp : substructure.getValue().met_positions) {
                    atom_states.get(sp.molecule_number).set(sp.atom_number, true);
                }
            }
        }
    }

    String drawMolecule(String inchi) throws Exception {

        // Generate factory - throws CDKException if native code does not load
        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
// Get InChIToStructure
        InChIToStructure intostruct = factory.getInChIToStructure(
                inchi, DefaultChemObjectBuilder.getInstance()
        );

        IAtomContainer container = intostruct.getAtomContainer();

        IMolecule molecule = new Molecule(container);

        SmilesGenerator sg = new SmilesGenerator();
        String name = sg.createSMILES(molecule);
        // layout the molecule
        StructureDiagramGenerator sdg = new StructureDiagramGenerator();
        sdg.setMolecule(molecule, false);
        sdg.generateCoordinates();

        // make generators
        List<IGenerator<IAtomContainer>> generators = new ArrayList<IGenerator<IAtomContainer>>();
        generators.add(new BasicSceneGenerator());
        generators.add(new BasicBondGenerator());
        //generators.add(new RingPlateGenerator());
        generators.add(new BasicAtomGenerator());

        // setup the renderer
        AtomContainerRenderer renderer = new AtomContainerRenderer(generators, new AWTFontManager());
        RendererModel model = renderer.getRenderer2DModel();
        model.set(BasicAtomGenerator.CompactAtom.class, true);
        model.set(BasicAtomGenerator.AtomRadius.class, 6.0);
        model.set(BasicAtomGenerator.CompactShape.class, BasicAtomGenerator.Shape.OVAL);
        model.set(BasicAtomGenerator.KekuleStructure.class, true);
        model.set(BasicBondGenerator.BondWidth.class, 5.0);

        // get the image
        Image image = new BufferedImage(400, 400, BufferedImage.TYPE_4BYTE_ABGR);
        Graphics2D g = (Graphics2D)image.getGraphics();
        g.setColor(Color.WHITE);
        g.fill(new Rectangle2D.Double(0, 0, 400, 400));

        // paint
        // renderer.paint(molecule, new AWTDrawVisitor(g));
        renderer.paint(molecule, new AWTDrawVisitor(g), new Rectangle2D.Double(0, 0, 400, 400), true);
        g.dispose();

        // write to file
        File file = new File("images", name + ".png");
        ImageIO.write((RenderedImage)image, "PNG", file);
        return name;
    }

    void dump_data(String path) throws Exception {
        PrintWriter writer = new PrintWriter(path + "final_signatures_active.txt", "UTF-8");
        for (Map.Entry<String, Double> substructure : final_signatures_active.entrySet()) {
            writer.println(substructure.getKey() + ' ' + substructure.getValue());
        }
        writer.close();
        writer = new PrintWriter(path + "final_signatures_inactive.txt", "UTF-8");
        for (Map.Entry<String, Double> substructure : final_signatures_inactive.entrySet()) {
            writer.println(substructure.getKey() + ' ' + substructure.getValue());
        }
        writer.close();
        writer = new PrintWriter(path + "signatures_active_heights.txt", "UTF-8");
        for (Map.Entry<String, Integer> substructure : heights_active.entrySet()) {
            writer.println(substructure.getKey() + ' ' + substructure.getValue());
        }
        writer.close();
        writer = new PrintWriter(path + "signatures_inactive_heights.txt", "UTF-8");
        for (Map.Entry<String, Integer> substructure : heights_inactive.entrySet()) {
            writer.println(substructure.getKey() + ' ' + substructure.getValue());
        }
        writer.close();
        writer = new PrintWriter("int_signatures_active.txt", "UTF-8");
        System.out.println(int_active.size());
        Iterator it = int_active.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            writer.println(drawMolecule((String)pair.getKey()));
            writer.println(pair.getValue().toString());
        }
        writer.close();
        writer = new PrintWriter("int_signatures_inactive.txt", "UTF-8");
        it = int_inactive.entrySet().iterator();
        System.out.println(int_inactive.size());
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            writer.println(drawMolecule((String)pair.getKey()));
            writer.println(pair.getValue().toString());
        }
        writer.close();
        System.out.println(final_signatures_active.size());
        System.out.println(final_signatures_inactive.size());
    }

    void run(String path) throws Exception {
        for (int height = 0; height < 20; height++) { //then redo while there are interesting substructures
            Map<String, Signature_state> current_signatures = new HashMap<String, Signature_state>();
            Map<String, String> signature_to_inchi = new HashMap<String, String>();
            find_substructures_with_height(current_signatures, signature_to_inchi, height);
            process_signatures(current_signatures, height);
        }
        dump_data(path);
    }

    double p_value_threshold;
    int minimum_occurrence;
    double substructure_frequency;
    int data_num;
    int active_num;
    double pi;
    boolean is_next;
    ArrayList<IMolecule> molecule_data;
    ArrayList<Boolean> activity_data;
    ArrayList<ArrayList<Boolean>> atom_states;
    HashMap<String, Double> final_signatures_active;
    HashMap<String, Double> final_signatures_inactive;
    HashMap<String, String> active;
    HashMap<String, String> inactive;
    HashMap<String, Double> prev_active;
    HashMap<String, Double> prev_inactive;
    HashMap<String, p_value_pair> int_active;
    HashMap<String, p_value_pair> int_inactive;
    HashMap<String, Integer> heights_active;
    HashMap<String, Integer> heights_inactive;

    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException, CDKException,
            IOException {
        //Algorithm a = new Algorithm(0.05, 5, 0.8);
        //a.data_initialization("pubchem_cyp1a2_train.sdf", "Inhibitor", "Noninhibitor");
        //a.run("");
        Observer observer = new Observer("ames_levenberg_experimental.sdf", "ames_levenberg_predicted.sdf");
        observer.initialize_algorithms(0.05, 5, 0.8);
        try {
            observer.run();
        } catch (Exception e) {
            e.printStackTrace();
        }
        observer.analyze();
        //Classification classifier;
        //classifier = new Classification();
        //classifier.data_initialization("");
        //classifier.classify_data("pubchem_cyp1a2_test.sdf");
        //a.classify_data("final_signatures_active.txt", "final_signatures_inactive.txt", "test_mol_for_class.sdf");
    }
}
