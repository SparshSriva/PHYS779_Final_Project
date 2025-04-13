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
                measurement_label = symbols("m" + str(measurement_counter))
                in1,in2,out1 = create_smooth_merge(graph_r, q_in1, q_in2, q_out1, measurement_label)
                measurement_counter += 1
                
                ## TODO: 
                #       -- Validate that output of merge is the same qubit as either input
                ##      -- Loop through arrays of inputs and outputs of self instead of num_qubits to set correct values
                
                for i in range(self.num_qubits):
                    if i != q_in1 and i != q_in2:
                        
                        max_row = graph_r.depth()
                        
                        in_qubit = graph_r.add_vertex(VertexType.BOUNDARY, qubit=i, row=0) 
                        out_qubit = graph_r.add_vertex(VertexType.BOUNDARY, qubit=i, row=max_row)
                        graph_r.add_edge((in_qubit, out_qubit))
                        temp_inputs.append(in_qubit)
                        temp_outputs.append(out_qubit)
                        
                temp_inputs.append(in1)
                temp_inputs.append(in2)
                temp_outputs.append(out1)
            
                graph_r.set_inputs((temp_inputs))
                graph_r.set_outputs((temp_outputs))
            
                print("Total graph: ",self.total_graph.outputs())
                print("graph_r: ",graph_r.inputs())
            
                self.total_graph.compose(graph_r)
                
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


def create_smooth_merge(graph: Graph, in1: int, in2: int, out: int, meas_outcome, pauli_frame: PauliFrame=None):
    """Create a smooth merge (Z-spider) with X correction on upper input leg.
    
    Args:
        graph: ZX graph to build upon
        in1: Qubit index for upper input boundary
        in2: Qubit index for lower input boundary
        out: Qubit index for output boundary
        meas_outcome: Symbolic measurement outcome
        pauli_frame: Optional PauliFrame object to record the X byproduct
    """

    inp1 = graph.add_vertex(VertexType.BOUNDARY, qubit=in1, row=0)
    inp2 = graph.add_vertex(VertexType.BOUNDARY, qubit=in2, row=0)
    correction = graph.add_vertex(VertexType.X, phase=meas_outcome*pi, qubit=in1, row=1)
    merge = graph.add_vertex(VertexType.Z, phase=0, qubit=in1, row=2)
    out_vert = graph.add_vertex(VertexType.BOUNDARY, qubit=out, row=3)

    graph.add_edges([(inp1, correction), (correction, merge),
                     (inp2, merge), (merge, out_vert)])

    graph.set_inputs([inp1, inp2])
    graph.set_outputs([out_vert])

    # If we want to track the correction in a PauliFrame:
    if pauli_frame is not None:
        pauli_frame.update(in1, "X")  # record X byproduct on qubit in1

    return inp1, inp2, out_vert


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