import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;

import java.io.*;

import java.math.BigInteger;

import java.util.*;
import java.util.List;

import com.google.common.math.BigIntegerMath;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecule;
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
              HashMap<String, Double> ex_active, HashMap<String, Double> ex_inactive, boolean original) {
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
        //this.heights_active = new HashMap<String, Integer>();
        //this.heights_inactive = new HashMap<String, Integer>();
        this.data_num = 0;
        this.active_num = 0;
        this.prev_active = ex_active;
        this.prev_inactive = ex_inactive;
        this.int_active = new HashMap<String, p_value_pair >();
        this.int_inactive = new HashMap<String, p_value_pair>();
        this.is_original = original;
        this.accuracies = new ArrayList<Double>();
    }

    void read_molecules(String pathname, String activity_name, String inactivity_name) {
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
            System.err.println("Original molecule size" + molecule_data.size());
            //here fill molecule_data, activity_data and allocate memory for atom states and fill it with 0
            //while initializing determine how many activity compounds there are
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
    }


    void read_molecules_predicted(String pathname, String activity_name, String inactivity_name) {
        try {
            File sdfFile = new File(pathname);
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(sdfFile), DefaultChemObjectBuilder.getInstance());
            int inactive_num = 0;
            int index = 0;
            while (reader.hasNext()) {
                IMolecule molecule = (IMolecule) reader.next();
                if (!activity_data.get(index)) {
                    ++inactive_num;
                    molecule_data.add(molecule);
                    ArrayList<Boolean> atom_states_local = new ArrayList<Boolean>(molecule.getAtomCount());
                    for (int t = 0; t < molecule.getAtomCount(); t++) {
                        atom_states_local.add(false);
                    }
                    atom_states.add(atom_states_local);
                    data_num++;
                } else {
                    ++active_num;
                    molecule_data.add(molecule);
                    ArrayList<Boolean> atom_states_local = new ArrayList<Boolean>(molecule.getAtomCount());
                    for (int t = 0; t < molecule.getAtomCount(); t++) {
                        atom_states_local.add(false);
                    }
                    atom_states.add(atom_states_local);
                    data_num++;
                }
                ++index;
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
        if (accuracies.size() != molecule_data.size()) {
            System.out.println("Non-matching sizes");
            System.out.println(accuracies.size());
            System.out.println(molecule_data.size());
        }
    }

    void read_accuracies(String pathname, String activity_name, String inactivity_name) throws Exception{
        Scanner scan;
        File file = new File(pathname.replace("_experimental.sdf", "_accuracies.csv"));
        try {
            scan = new Scanner(file);
            scan.nextLine();

            while(scan.hasNextLine()) {
                String[] data = scan.nextLine().split(";");
                if (data.length == 0) {
                    activity_data.add(false);
                    this.accuracies.add(0.0);
                    continue;
                }
                if (!is_original) {
                    if (data[0].equals(activity_name)) {
                        activity_data.add(true);
                    } else {
                        activity_data.add(false);
                    }
                }
                this.accuracies.add(Double.parseDouble(data[2].replaceAll(",", ".")));
            }
            System.out.println(activity_data.size());

        } catch (FileNotFoundException e1) {
            e1.printStackTrace();
        }
        /*if (accuracies.size() != molecule_data.size()) {
            throw new Exception("Unmatching accuracies and data file");
        }*/
    }

    void data_initialization(String pathname, String activity_name, String inactivity_name, boolean use_accuracies) throws Exception{
        use_accuracy = use_accuracies;
        if (!is_original && use_accuracies) {
            read_accuracies(pathname, activity_name, inactivity_name);
            read_molecules_predicted(pathname, activity_name, inactivity_name);
        } else {
            read_molecules(pathname, activity_name, inactivity_name);
        }
    }

    class Signature_position {
        int molecule_number;
        int atom_number;

        Signature_position(int i, int j) {
            this.molecule_number = i;
            this.atom_number = j;
        }
    }

    class Signature_state {
        Signature_state() {
            this.met_positions = new ArrayList<Signature_position>();
            n_dash = 0;
            m_dash = 0;
            average_accuracy = 0;
        }

        void add(int i, int j, Boolean activity, double accuracy) {
            Signature_position sp = new Signature_position(i, j);
            this.met_positions.add(sp);
            n_dash++;
            if (activity)
                m_dash++;
            average_accuracy += accuracy;
        }

        //int processed_condition; //not met; met but not interesting; met but interesting; included in final list
        ArrayList<Signature_position> met_positions; //positions of atoms, where this signature was met
        double n_dash;//number of substructures
        double m_dash;//number of active substructures
        double average_accuracy; //average accuracy of predicted activity
    }

    class SignatureInfo {
        double active_entry_count_original;
        double inactive_entry_count_original;
        double active_entry_count_trained;
        double inactive_entry_count_trained;
        double p_value_original;
        double p_value_trained;
        double average_accuracy;
        int category;
        String smiles;
        SignatureInfo() {
            active_entry_count_original = 0;
            inactive_entry_count_original = 0;
            active_entry_count_trained = 0;
            inactive_entry_count_trained = 0;
            p_value_original = 0;
            p_value_trained = 0;
            category = 10;
        }
        void fill_original(double active, double inactive, double p_value) {
            active_entry_count_original = active;
            inactive_entry_count_original = inactive;
            p_value_original = p_value;
        }
        void fill_trained(double active, double inactive, double p_value, double accuracy) {
            active_entry_count_trained = active;
            inactive_entry_count_trained = inactive;
            p_value_trained = p_value;
            average_accuracy = accuracy / (active + inactive);
        }
        void fill_category(int input_category) {
            category = input_category;
        }
        void fill_smiles(String smiles) {
            this.smiles = smiles;
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
            HashMap<String, HashSet<Integer>> inchi_indices = new HashMap<String, HashSet<Integer>>();
            for (int i = 0; i < m.getAtomCount(); i++) {
                HashSet<Integer> local_hash_set = new HashSet<Integer>();
                if (this.atom_states.get(molecule_index).get(i) == false) { //not in the current signatures and not viewed (we still search in this direction and it's not viewed)
                    AtomSignature as = new AtomSignature(i, height, m);
                    String string_repr = atom_signature_to_inchi(as);
                    int vertex_count = as.getVertexCount();
                    for (int iii = 0; iii < vertex_count; ++iii) {
                        local_hash_set.add(as.getOriginalVertexIndex(iii));
                    }
                    if (inchi_indices.get(string_repr) == null) {
                        inchi_indices.put(string_repr, local_hash_set);
                    }
                    else {
                        if (inchi_indices.get(string_repr).containsAll(local_hash_set)) {
                            Signature_position sp = new Signature_position(molecule_index, i);
                            current_signatures.get(string_repr).met_positions.add(sp);
                            continue;
                        } else {
                            inchi_indices.get(string_repr).addAll(local_hash_set);
                        }
                    }
                    if (current_signatures.containsKey(string_repr) == false) {
                        Signature_state ss = new Signature_state();
                        current_signatures.put(string_repr, ss);
                    }
                    if (!is_original) {
                        if (use_accuracy) {
                            current_signatures.get(string_repr).add(molecule_index, i,
                                                                    this.activity_data.get(molecule_index),
                                                                    accuracies.get(molecule_index));
                        }
                        else {
                            current_signatures.get(string_repr).add(molecule_index, i,
                                                                    this.activity_data.get(molecule_index), 0);
                        }
                    } else {
                        current_signatures.get(string_repr).add(molecule_index, i,
                                                                this.activity_data.get(molecule_index), 0);
                    }
                }
            }
            molecule_index++;
        }
    }

    Double calculate_p_value(double m_, double n_) {
        Double true_p_value = 0.0;
        for (double k = m_; k < n_; k++) {
            if (n_ - k == 0 || k == 0 || data_num == 0) {
                System.out.println(n_ - k + ' ' + k + ' ' + data_num);
            }
            BigInteger coef = BigIntegerMath.binomial((int) n_, (int) k);
            true_p_value += coef.doubleValue() * Math.pow(new Double(active_num)/data_num, k) *
                    Math.pow(new Double(data_num - active_num)/data_num, n_ - k);
        }
        true_p_value += Math.pow(new Double(active_num) / data_num, n_);
        return true_p_value;
    }

    void process_signatures(Map<String, Signature_state> current_signatures, int height, Map<String, SignatureInfo> table)
            throws CDKException, FileNotFoundException, UnsupportedEncodingException {
        for (Map.Entry<String, Signature_state> substructure : current_signatures.entrySet()) {
            String as_string = substructure.getKey();
            if (substructure.getValue().n_dash > minimum_occurrence) {
                double current_frequency = substructure.getValue().m_dash / substructure.getValue().n_dash;
                double n_ = substructure.getValue().n_dash;
                double m_ = substructure.getValue().m_dash;
                if (current_frequency > substructure_frequency) {
                    Double true_p_value = calculate_p_value(m_, n_);
                    if (true_p_value < p_value_threshold) {
                        for (Signature_position sp : substructure.getValue().met_positions) {
                            atom_states.get(sp.molecule_number).set(sp.atom_number, true);
                        }
                        final_signatures_active.put(as_string, true_p_value);
                        active.put(as_string, as_string);
                        if (is_original) {
                            SignatureInfo si = new SignatureInfo();
                            si.fill_original(m_, n_ - m_, true_p_value);
                            table.put(as_string, si);
                        }
                        else {
                            SignatureInfo si = table.get(as_string);
                            if (si == null) {
                                si = new SignatureInfo();
                                table.put(as_string, si);
                            }
                            si.fill_trained(m_, n_ - m_, true_p_value, substructure.getValue().average_accuracy);
                        }
                    }
                    else {
                        if (prev_active.size() != 0) {
                            if (prev_active.containsKey(as_string)) {
                                int_active.put(as_string, new p_value_pair(prev_active.get(as_string), true_p_value));
                            }
                        }
                    }
                }

                current_frequency = (substructure.getValue().n_dash - substructure.getValue().m_dash) / substructure.getValue().n_dash;

                if (current_frequency > substructure_frequency) {
                    Double true_p_value = 0.0;
                    m_ = n_ - substructure.getValue().m_dash;
                    for (double k = m_; k < n_; k++) {
                        BigInteger coef = BigIntegerMath.binomial((int) n_, (int) k);
                        true_p_value += coef.doubleValue() * Math.pow(new Double(active_num)/data_num, k) *
                                Math.pow(new Double(data_num - active_num)/data_num, n_ - k);
                    }
                    true_p_value += Math.pow(new Double(active_num) / data_num, n_);
                    if (as_string.equals("[C]([C]([C])[C]([C][C]))")) {
                        System.out.println("Found interesting " + substructure.getValue().m_dash + " " + substructure.getValue().n_dash + " " + height + " " + true_p_value);
                    }
                    if (true_p_value < p_value_threshold) {
                        for (Signature_position sp : substructure.getValue().met_positions) {
                            atom_states.get(sp.molecule_number).set(sp.atom_number, true);
                        }
                        final_signatures_inactive.put(as_string, true_p_value);
                        inactive.put(as_string, as_string);
                        if (is_original) {
                            SignatureInfo si = new SignatureInfo();
                            si.fill_original(substructure.getValue().m_dash, n_ - substructure.getValue().m_dash, true_p_value);
                            table.put(as_string, si);
                        }
                        else {
                            SignatureInfo si = table.get(as_string);
                            if (si == null) {
                                si = new SignatureInfo();
                                table.put(as_string, si);
                            }
                            si.fill_trained(substructure.getValue().m_dash, n_ - substructure.getValue().m_dash, true_p_value, substructure.getValue().average_accuracy);
                        }
                    }
                    else {
                        if (prev_inactive.size() != 0) {
                            if (prev_inactive.containsKey(as_string)) {
                                int_inactive.put(as_string, new p_value_pair(prev_inactive.get(as_string), true_p_value));
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
        try {
            sdg.generateCoordinates();
        }
        catch (Exception e) {
            System.err.println(inchi);
            return "";
        }

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
        name = name.replace("#", "?");
        //System.out.println(name);
        File file = new File("images", name + ".png");
        ImageIO.write((RenderedImage)image, "PNG", file);
        return name;
    }

    void dump_data(String path) throws Exception {
        PrintWriter writer = new PrintWriter(path + "final_signatures_active.txt", "UTF-8");
        for (Map.Entry<String, Double> substructure : final_signatures_active.entrySet()) {
            try {
                if (drawMolecule(substructure.getKey().toString()) != "")
                    writer.println(substructure.getKey() + ' ' + substructure.getValue());
            } catch (Exception e) {
                if (substructure.getKey() != null) {
                    System.err.println(substructure.getKey().toString());
                } else {
                    System.err.println(e.getMessage());
                }
            }
        }
        writer.close();
        writer = new PrintWriter(path + "final_signatures_inactive.txt", "UTF-8");
        for (Map.Entry<String, Double> substructure : final_signatures_inactive.entrySet()) {
            try {
                if (drawMolecule(substructure.getKey().toString()) != "")
                    writer.println(substructure.getKey() + ' ' + substructure.getValue());
            } catch (Exception e) {
                if (substructure.getKey() != null) {
                    System.err.println(substructure.getKey().toString());
                } else {
                    System.err.println(e.getMessage());
                }
            }
        }
        PrintWriter writer2 = new PrintWriter("int_signatures_active.txt", "UTF-8");
        Iterator it = int_active.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            writer2.println(pair.getValue().toString());
        }
        writer2.close();
        writer2 = new PrintWriter("int_signatures_inactive.txt", "UTF-8");
        it = int_inactive.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            writer2.println(pair.getValue().toString());
        }
        writer2.close();
        System.out.println("Final signatures " + final_signatures_active.size());
        System.out.println("Final signatures " + final_signatures_inactive.size());
    }

    void run(Map<String, SignatureInfo> table, String path) throws Exception {
        for (int height = 0; height < 20; height++) { //then redo while there are interesting substructures
            Map<String, Signature_state> current_signatures = new HashMap<String, Signature_state>();
            Map<String, String> signature_to_inchi = new HashMap<String, String>();
            find_substructures_with_height(current_signatures, signature_to_inchi, height);
            process_signatures(current_signatures, height, table);
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
    boolean is_original;
    boolean use_accuracy;
    ArrayList<IMolecule> molecule_data;
    ArrayList<Boolean> activity_data;
    ArrayList<ArrayList<Boolean>> atom_states;
    ArrayList<Double> accuracies; //accuracies for predicted values
    HashMap<String, Double> final_signatures_active;
    HashMap<String, Double> final_signatures_inactive;
    HashMap<String, String> active;
    HashMap<String, String> inactive;
    HashMap<String, Double> prev_active;
    HashMap<String, Double> prev_inactive;
    HashMap<String, p_value_pair> int_active;
    HashMap<String, p_value_pair> int_inactive;

    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException, CDKException,
            IOException {
        Observer observer = new Observer(args[0], args[1], args[2], args[3]);
        Map<String, SignatureInfo> table = new HashMap<String, SignatureInfo>();

        try {
            observer.initialize_algorithms(Double.parseDouble(args[4]), Integer.parseInt(args[5]), Double.parseDouble(args[6]), Boolean.parseBoolean(args[7]));
            observer.run(table);
        } catch (Exception e) {
            e.printStackTrace();
        }
        observer.analyze(table);
    }
}
