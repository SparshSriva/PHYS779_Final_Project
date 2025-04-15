# ZXLS.py

from pyzx import Graph, VertexType, simplify, draw
from sympy import symbols, pi


class PauliFrame:
    """Tracks the classical record of byproduct corrections for each qubit."""

    def __init__(self):
        self.frame = {}  # Maps qubit index -> Pauli ("I", "X", "Z", or "Y")

    def update(self, qubit: int, new: str):
        """Compose a new Pauli with the existing one for 'qubit'."""
        old = self.frame.get(qubit, "I")
        self.frame[qubit] = self._compose_paulis(old, new)

    def _compose_paulis(self, a: str, b: str) -> str:
        """Table for multiplying two Pauli operators ignoring phases."""
        pauli_table = {
            ("I", "I"): "I", ("I", "X"): "X", ("I", "Z"): "Z", ("I", "Y"): "Y",
            ("X", "I"): "X", ("X", "X"): "I", ("X", "Z"): "Y", ("X", "Y"): "Z",
            ("Z", "I"): "Z", ("Z", "X"): "Y", ("Z", "Z"): "I", ("Z", "Y"): "X",
            ("Y", "I"): "Y", ("Y", "X"): "Z", ("Y", "Z"): "X", ("Y", "Y"): "I"
        }
        return pauli_table[(a, b)]

    def __repr__(self):
        return f"PauliFrame({self.frame})"

class ParseInstructions:
    def __init__(self, num_qubits):
        if not isinstance(num_qubits, int) or num_qubits <= 0:
            raise ValueError("Number of qubits must be a positive integer")

        operation_type_lookup = {
            "s-split": {"type": "split", "mode": "s"},
            "r-split": {"type": "split", "mode": "r"},
            "s-merge": {"type": "merge", "mode": "s"},
            "r-merge": {"type": "merge", "mode": "r"},
            "x-measure": {"type": "measurement", "mode":"x-basis"},
            "z-measure": {"type": "measurement", "mode":"z-basis"},
            "init": {"type": "initialize", "mode": "ancilla"}
        }
        
        self.num_qubits = num_qubits
        self.total_graph = Graph()
        self.operations = []
        self.operation_type = operation_type_lookup
        
        self.inputs = []
        self.outputs = []
        
        for i in range(self.num_qubits):
            in_qubit = self.total_graph.add_vertex(VertexType.BOUNDARY, qubit=i, row=0) # Control Input
            out_qubit = self.total_graph.add_vertex(VertexType.BOUNDARY, qubit=i, row=1)
            self.total_graph.add_edge((in_qubit, out_qubit))
            self.inputs.append(in_qubit)
            self.outputs.append(out_qubit)
            
        self.total_graph.set_inputs(self.inputs)
        self.total_graph.set_outputs(self.outputs)

    def validate_operation(self, operation):
        if operation in self.operation_type:
            return True
        else:
            return False

    def add_operation(self, operation, *qubits):
        """Add operation with qubit indices"""
        if not self.validate_operation(operation):
            raise TypeError("invalid-operation")
        if not all(isinstance(q, int) and 0 <= q < self.num_qubits for q in qubits):
            raise ValueError("Qubit indices must be valid integers within the range")
        self.operations.append((operation, qubits))

    def get_non_matching_qubit_keys(self, g_qubits_dict, check_inputs=True):
        total_qubits = self.total_graph.qubits()

        if check_inputs:
            relevant_keys = self.total_graph.inputs()
        else:
            relevant_keys = self.total_graph.outputs()

        boundary_values = set(total_qubits[k] for k in relevant_keys if k in total_qubits)

        g_values = set(g_qubits_dict.values())

        missing_values = boundary_values - g_values
    
        uncovered_indices = [
            index for index in relevant_keys
            if total_qubits.get(index) in missing_values
        ]

        return uncovered_indices


    def execute_operations(self):
        
        measurement_counter = 0

        for operation, qubits in self.operations:
            graph_r = Graph()
            if operation == "s-split":
                if len(qubits) != 3:
                    raise ValueError("s-split operation requires 3 qubit indices: input, output1, output2")
                q_in, q_out1, q_out2 = qubits
                measurement_label = symbols("m" + str(measurement_counter))
                create_smooth_split(graph_r, q_in, q_out1, q_out2, measurement_label)
                measurement_counter += 1
                #self.total_graph.compose(graph_r)
                
            elif operation == "s-merge":
                
                temp_inputs = []
                temp_outputs = []
            
            
                if len(qubits) != 3:
                    raise ValueError("s-merge operation requires 3 qubit indices: input1, input2, output1")
                q_in1, q_in2, q_out1 = qubits
                
                if q_out1 != q_in1 and q_out1 != q_in2:
                    raise ValueError("outputs must be equal to one the inputs")
                
                measurement_label = symbols("m" + str(measurement_counter))
                in1,in2,out1 = create_smooth_merge(graph_r, self.total_graph, q_in1, q_in2, q_out1, measurement_label)
                measurement_counter += 1
                
                ## Maybe fix?? --> Loop through arrays of inputs and outputs of self instead of num_qubits to set correct values
                
                temp_inputs.append(in1)
                temp_inputs.append(in2)
                temp_outputs.append(out1)
                
                print("total g inputs: ", self.total_graph.inputs())
                print("total g outputs: ", self.total_graph.outputs())
                print("total g qubits: ", self.total_graph.qubits())
                
                print("g inputs: ", graph_r.inputs())
                print("g outputs: ", graph_r.outputs())
                print("g qubits: ", graph_r.qubits())
                
                non_matching = self.get_non_matching_qubit_keys(graph_r.qubits(), True)
                
                for i in non_matching:
                    if i not in graph_r.vertex_set():
                        graph_r.add_vertex_indexed(i)
                        graph_r.set_type(i, VertexType.BOUNDARY)
                        graph_r.set_row(i, 0)
                        graph_r.set_qubit(i, self.total_graph.qubit(i))
                        temp_inputs.append(i)
                    
                ## Fix this STOOPID, add all the outputs that dont have a respective input from the subgraph
                non_matching = self.get_non_matching_qubit_keys(graph_r.qubits(), False)
                for i in non_matching:
                    if i not in graph_r.vertex_set():
                        graph_r.add_vertex_indexed(i)
                        graph_r.set_type(i, VertexType.BOUNDARY)
                        graph_r.set_row(i, graph_r.depth())
                        graph_r.set_qubit(i, self.total_graph.qubit(i))
                        temp_outputs.append(i)
                
                draw(graph_r, labels=True)
                
                for i in self.total_graph.inputs():
                    qubit_validate = graph_r.qubit(i)
                    if i not in graph_r.inputs():
                        graph_r.add_edge((i,i+1))
            
                graph_r.set_inputs((temp_inputs))
                graph_r.set_outputs((temp_outputs))
            
                self.total_graph.compose(graph_r)
            
            elif operation == "r-merge":
                
                temp_inputs = []
                temp_outputs = []
            
            
                if len(qubits) != 3:
                    raise ValueError("r-merge operation requires 3 qubit indices: input1, input2, output1")
                q_in1, q_in2, q_out1 = qubits
                
                if q_out1 != q_in1 and q_out1 != q_in2:
                    raise ValueError("outputs must be equal to one the inputs")
                
                measurement_label = symbols("m" + str(measurement_counter))
                in1,in2,out1 = create_rough_merge(graph_r, q_in1, q_in2, q_out1, measurement_label)
                measurement_counter += 1
                
                ## Maybe fix?? --> Loop through arrays of inputs and outputs of self instead of num_qubits to set correct values
                
                temp_inputs.append(in1)
                temp_inputs.append(in2)
                

                
                for input in self.total_graph.inputs():
                    for output in self.total_graph.outputs():
                        if input != q_in1 and input != q_in2:
                            temp_inputs.append(input)
                            temp_outputs.append(output)
                        
                temp_outputs.append(out1)
            
                graph_r.set_inputs((temp_inputs))
                graph_r.set_outputs((temp_outputs))
            
                self.total_graph.compose(graph_r)   
                
            else:
                raise ValueError("invalid operation")
                
    def list_operations(self):
        return self.operations

    def __repr__(self):
        return f"ParseInstructions(num_qubits={self.num_qubits}, operations={self.operations})"

def create_smooth_split(graph: Graph, inp: int, out1: int, out2: int, meas_outcome):
    """Create a smooth split (Z-spider) with X correction on upper output leg.
    
    Args:
        graph: ZX-Graph to build upon
        inp: Qubit index for input boundary
        out1: Qubit index for upper output boundary
        out2: Qubit index for lower output boundary
        meas_outcome: Symbolic measurement outcome (e.g., m1)
        
    Returns:
        Tuple of (upper_boundary_vertex, lower_boundary_vertex, input_boundary_vertex)
    """

    inp_vert = graph.add_vertex(VertexType.BOUNDARY, qubit=inp, row=0)
    split = graph.add_vertex(VertexType.Z, phase=0, qubit=inp, row=1)
    correction = graph.add_vertex(VertexType.X, phase=meas_outcome*pi, qubit=inp, row=2)
    upper = graph.add_vertex(VertexType.BOUNDARY, qubit=out1, row=3)
    lower = graph.add_vertex(VertexType.BOUNDARY, qubit=out2, row=3)

    # Connect edges
    graph.add_edges([(inp_vert, split), (split, correction),
                     (correction, upper), (split, lower)])

    graph.set_inputs([inp_vert])
    graph.set_outputs([upper, lower])
    return upper, lower, inp_vert


def create_rough_split(graph: Graph, inp: int, out1: int, out2: int, meas_outcome):
    """Create a rough split (X-spider) with Z correction on upper output leg.
    
    Args:
        graph: ZX-Graph to build upon
        inp: Qubit index for input boundary
        out1: Qubit index for upper output boundary
        out2: Qubit index for lower output boundary
        meas_outcome: Symbolic measurement outcome (e.g., m1)
        
    Returns:
        Tuple of (upper_boundary_vertex, lower_boundary_vertex, input_boundary_vertex)
    """

    inp_vert = graph.add_vertex(VertexType.BOUNDARY, qubit=inp, row=0)
    split = graph.add_vertex(VertexType.X, phase=0, qubit=inp, row=1)
    correction = graph.add_vertex(VertexType.Z, phase=meas_outcome*pi, qubit=out1, row=2)
    upper = graph.add_vertex(VertexType.BOUNDARY, qubit=out1, row=3)
    lower = graph.add_vertex(VertexType.BOUNDARY, qubit=out2, row=3)

    graph.add_edges([(inp_vert, split), (split, correction),
                     (correction, upper), (split, lower)])

    graph.set_inputs([inp_vert])
    graph.set_outputs([upper, lower])
    return upper, lower, inp_vert


def create_smooth_merge(graph: Graph, full_graph: Graph, in1: int, in2: int, out: int, meas_outcome, pauli_frame: PauliFrame=None):
    """Create a smooth merge (Z-spider) with X correction on upper input leg.
    
    Args:
        graph: ZX graph to build upon
        in1: Qubit index for upper input boundary
        in2: Qubit index for lower input boundary
        out: Qubit index for output boundary
        meas_outcome: Symbolic measurement outcome
        pauli_frame: Optional PauliFrame object to record the X byproduct
    """
    
    all_qubits = full_graph.qubits()
    
    for index in full_graph.outputs():
        if index in all_qubits.keys():
            if full_graph.qubit(index) == in1:
                inp1 = index
            if full_graph.qubit(index) == in2:
                inp2 = index
            if full_graph.qubit(index) == out:
                outp = index


    q_in_1 = full_graph.qubit(inp1)
    q_in_2 = full_graph.qubit(inp2)

    ## Input 1
    graph.add_vertex_indexed(inp1)
    graph.set_type(inp1, VertexType.BOUNDARY)
    ## TODO: Figure out logic to set row
    ## Maybe this works? Test further
    graph.set_row(inp1, full_graph.row(q_in_1))
    ## TODO: Set index AND qubit number here
    ## Maybe this works? Test further
    graph.set_qubit(inp1, q_in_1)
    
    ## Output
    n_out = graph.add_vertex(VertexType.BOUNDARY, qubit=out, row=full_graph.row(q_in_1)+3)

    ## Input 2
    graph.add_vertex_indexed(inp2)
    graph.set_type(inp2, VertexType.BOUNDARY)
    ## TODO: Figure out logic to set row
    ## Maybe this works? Test further
    graph.set_row(inp2, full_graph.row(q_in_1))
    ## TODO: Set index AND qubit number here
    ## Maybe this works? Test further
    graph.set_qubit(inp2, q_in_2)
    
    total_vertex_count = len(full_graph.vertex_set())
    
    ## Correction  correction = graph.add_vertex(VertexType.X, phase=meas_outcome*pi, qubit=in1, row=inp_row+1)
    graph.add_vertex_indexed(total_vertex_count)
    graph.set_type(total_vertex_count, VertexType.X)
    graph.set_phase(total_vertex_count, meas_outcome*pi)
    graph.set_row(total_vertex_count, full_graph.row(q_in_1)+1)
    graph.set_qubit(total_vertex_count, in1)
    
    ## Merge merge = graph.add_vertex(VertexType.Z, phase=0, qubit=in1, row=inp_row+2)
    graph.add_vertex_indexed(total_vertex_count+1)
    graph.set_type(total_vertex_count+1, VertexType.Z)
    graph.set_phase(total_vertex_count+1, 0)
    graph.set_row(total_vertex_count+1, full_graph.row(q_in_1)+2)
    graph.set_qubit(total_vertex_count+1, in1)
    

    graph.add_edges([(inp1, total_vertex_count), (total_vertex_count, total_vertex_count+1),
                     (inp2, total_vertex_count+1), (total_vertex_count+1, n_out)])

    graph.set_inputs([in1, in2])
    graph.set_outputs([n_out])

    # If we want to track the correction in a PauliFrame:
    if pauli_frame is not None:
        pauli_frame.update(in1, "X")  # record X byproduct on qubit in1

    return in1, in2, n_out


def create_rough_merge(graph: Graph, in1: int, in2: int, out: int, meas_outcome):
    """Create a rough merge (X-spider) with Z correction on upper input leg.
    
    Args:
        graph: ZX graph to build upon
        in1: Qubit index for upper input boundary
        in2: Qubit index for lower input boundary
        out: Qubit index for output boundary
        meas_outcome: Symbolic measurement outcome
    """

    inp1 = graph.add_vertex(VertexType.BOUNDARY, qubit=in1, row=0)
    inp2 = graph.add_vertex(VertexType.BOUNDARY, qubit=in2, row=0)
    correction = graph.add_vertex(VertexType.Z, phase=meas_outcome*pi, qubit=in1, row=1)
    merge = graph.add_vertex(VertexType.X, phase=0, qubit=in1, row=2)
    out_vert = graph.add_vertex(VertexType.BOUNDARY, qubit=out, row=3)

    graph.add_edges([(inp1, correction), (correction, merge),
                     (inp2, merge), (merge, out_vert)])

    graph.set_inputs([inp1, inp2])
    graph.set_outputs([out_vert])
    return inp1, inp2, out_vert


def substitute_measurements(graph: Graph, substitution_dict: dict):
    """Substitute symbolic measurement results into the graph and simplify it.
    
    For each vertex in the graph, if its phase is a Sympy expression, we substitute the
    values in substitution_dict; otherwise, we leave it unchanged.
    
    Args:
        graph: A ZX Graph (e.g., of type GraphS) containing symbolic phases.
        substitution_dict: A dictionary mapping symbols (e.g., m1) to numerical values (e.g., 0 or 1).
        
    Returns:
        A simplified copy of the graph with the substitutions applied.
    """
    from sympy import Basic
    g_sub = graph.copy()
    for v in list(g_sub.vertices()):
        ph = g_sub.phase(v)
        if ph is not None:
            if isinstance(ph, Basic):
                new_ph = ph.subs(substitution_dict)
            else:
                new_ph = ph  # Already a number, no substitution needed
            g_sub.set_phase(v, new_ph)
    simplify.full_reduce(g_sub)
    return g_sub