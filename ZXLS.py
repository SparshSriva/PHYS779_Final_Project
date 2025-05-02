# ZXLS.py

from pyzx import Graph, VertexType, simplify, draw
from sympy import symbols, pi, Basic


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

def get_boundary_vertex_by_qubit(graph, boundary_list, qubit_value):
    """
    Returns the vertex ID from `boundary_list` that corresponds to `qubit_value`.
    """
    for v in boundary_list:
        if graph.qubit(v) == qubit_value:
            return v
    raise ValueError(f"No vertex in boundary list for qubit {qubit_value}")

def clean_orphan_boundaries(graph, temp_inputs, temp_outputs):
    for v in temp_inputs:
        if len(graph.neighbors(v)) == 0:
            graph.remove_vertex(v)
    for v in temp_outputs:
        if len(graph.neighbors(v)) == 0:
            graph.remove_vertex(v)
            
class ParseInstructions:
    def __get_qubits_not_involved(self, used_qubits):
        """
        Returns the dictionary of inputs and outputs (with their index) from the graph that do not appear in the used_qubits set.
        """
        qubits_inputs = {}
        qubits_outputs = {}
        for v in self.total_graph.inputs():
            if self.total_graph.qubit(v) not in used_qubits:
                qubits_inputs[v] = self.total_graph.qubit(v)

        for v in self.total_graph.outputs():
            if self.total_graph.qubit(v) not in used_qubits:
                qubits_outputs[v] = self.total_graph.qubit(v)

        return qubits_inputs, qubits_outputs
    
    def __init__(self, num_qubits,ancilla_qubits=None):
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
        ## TODO: verify ancilla qubits should be checked or not in other protocols
        self.ancilla_info = ancilla_qubits or []
        ancilla_lookup = {entry["qubit_number"]: entry for entry in self.ancilla_info}

        self.inputs = []
        self.outputs = []
        
        if len(self.ancilla_info) != 0:
            for i in self.ancilla_info:
                self._validate_ancilla_entry(i)
        
        for i in range(self.num_qubits):
            if i in ancilla_lookup:
                anc = ancilla_lookup[i]
                v = self.total_graph.add_vertex(
                    VertexType.Z if anc["basis"] == "Z" else VertexType.X,
                    qubit=i,
                    row=0,
                    phase=anc["phase"]
                )
                out_qubit = self.total_graph.add_vertex(VertexType.BOUNDARY, qubit=i, row=1)
                self.total_graph.add_edge((v,out_qubit))
                self.inputs.append(v)
                self.outputs.append(out_qubit)
            else:
                in_qubit = self.total_graph.add_vertex(VertexType.BOUNDARY, qubit=i, row=0)
                out_qubit = self.total_graph.add_vertex(VertexType.BOUNDARY, qubit=i, row=1)
                self.total_graph.add_edge((in_qubit, out_qubit))
                self.inputs.append(in_qubit)
                self.outputs.append(out_qubit)
                self.total_graph.set_inputs(self.inputs)
                self.total_graph.set_outputs(self.outputs)
    
    def _validate_ancilla_entry(self, entry):
        required_keys = {"qubit_number", "basis", "phase"}
        if not all(key in entry for key in required_keys):
            raise ValueError(f"Missing keys in ancilla definition: {entry}")
        if not isinstance(entry["qubit_number"], int):
            raise TypeError(f"Qubit number must be int, got {entry['qubit_number']}")
        if entry["basis"] not in ("X", "Z"):
            raise ValueError(f"Unsupported basis '{entry['basis']}' in ancilla definition")
        if entry["phase"] != 0 and entry["phase"] != 1:
            raise TypeError(f"Phase must be int 0 or 1, got {entry['phase']}")

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

    def get_uncovered_output_indices(self, subgraph):
        total_output_indices = set(self.total_graph.outputs())
        subgraph_output_indices = set(subgraph.outputs())
        return list(total_output_indices - subgraph_output_indices)

    def get_uncovered_input_indices(self, subgraph):
        total_output_indices = set(self.total_graph.inputs())
        subgraph_output_indices = set(subgraph.outputs())
        return list(total_output_indices - subgraph_output_indices)


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
                temp_inputs = []
                temp_outputs = []

                if len(qubits) != 3:
                    raise ValueError("s-split operation requires 3 qubit indices: input, output1, output2")
                q_in, q_out1, q_out2 = qubits

                if q_out1 != q_in:
                    raise ValueError("output 1 must be equal to the input")

                measurement_label = symbols("m" + str(measurement_counter))
                measurement_counter += 1

                upper, lower, inp_v = create_smooth_split(graph_r, self.total_graph, q_in, q_out1, q_out2, measurement_label)

                temp_inputs.append(inp_v)
                temp_outputs.extend([upper, lower])

                used_qubits = {q_out1, q_out2}
                
                unused_inputs, unused_outputs = self.__get_qubits_not_involved(used_qubits)

                for i in unused_inputs:
                    q = self.total_graph.qubit(i)
                    if q in used_qubits:
                        continue
                    if i not in graph_r.vertex_set():
                        graph_r.add_vertex_indexed(i)
                        graph_r.set_type(i, VertexType.BOUNDARY)
                        graph_r.set_row(i, 0)
                        graph_r.set_qubit(i, q)
                        temp_inputs.append(i)
                
                for i in unused_outputs:
                    q = self.total_graph.qubit(i)
                    if q in used_qubits:
                        continue
                    if i not in graph_r.vertex_set():
                        graph_r.add_vertex_indexed(i)
                        graph_r.set_type(i, VertexType.BOUNDARY)
                        graph_r.set_row(i, graph_r.depth())
                        graph_r.set_qubit(i, q)
                        temp_outputs.append(i)
                    
                temp_inputs.sort()
                temp_outputs.sort()
                
                graph_r.set_inputs(temp_inputs)
                graph_r.set_outputs(temp_outputs)    
                        
                for i in graph_r.outputs():
                    q = graph_r.qubit(i)
                    if q in used_qubits:
                        continue 
                    for j in graph_r.inputs():
                        if graph_r.qubit(j) == q:
                            if not graph_r.connected(j, i):
                                graph_r.add_edge((j, i))
                            break
                
                if len(graph_r.inputs()) != len(self.total_graph.outputs()):
                    print("graph_r inputs: ", len(graph_r.inputs()))
                    print("total_graph outputs: ", len(self.total_graph.outputs()))
                    raise ValueError("graph_r inputs do not match total_graph outputs")
                
                self.total_graph.compose(graph_r)

                
            elif operation == "init":
                continue
                
            elif operation == "s-merge":
                
                temp_inputs = []
                temp_outputs = []
            
            
                if len(qubits) != 3:
                    raise ValueError("s-merge operation requires 3 qubit indices: input1, input2, output1")
                q_in1, q_in2, q_out1 = qubits
                
                if q_out1 != q_in1 and q_out1 != q_in2:
                    raise ValueError("outputs must be equal to one of the inputs")
                
                measurement_label = symbols("m" + str(measurement_counter))
                in1,in2,out1 = create_smooth_merge(graph_r, self.total_graph, q_in1, q_in2, q_out1, measurement_label)
                measurement_counter += 1
                
                temp_inputs.append(in1)
                temp_inputs.append(in2)
                temp_outputs.append(out1)
                
                non_matching = self.get_non_matching_qubit_keys(graph_r.qubits(), True)
                
                used_qubits = {q_in1, q_in2}  

                for i in non_matching:
                    q = self.total_graph.qubit(i)
                    if q in used_qubits:
                        continue  
                    
                    if i not in graph_r.vertex_set():                        
                        graph_r.add_vertex_indexed(i)
                        graph_r.set_type(i, VertexType.BOUNDARY)
                        graph_r.set_row(i, 0)
                        graph_r.set_qubit(i, q)
                        temp_inputs.append(i)

                
                uncovered_output_indices = self.get_uncovered_output_indices(graph_r)

                for i in uncovered_output_indices:
                    q = self.total_graph.qubit(i)
                    if q in used_qubits:
                        continue  
                    
                    if i not in graph_r.vertex_set():
                        graph_r.add_vertex_indexed(i)
                        graph_r.set_type(i, VertexType.BOUNDARY)
                        graph_r.set_row(i, graph_r.depth())
                        graph_r.set_qubit(i, q)
                        temp_outputs.append(i)


                temp_inputs.sort()
                temp_outputs.sort()
                
                graph_r.set_inputs((temp_inputs))
                graph_r.set_outputs((temp_outputs))
                
                for i in graph_r.inputs():
                    q = graph_r.qubit(i)
                
                    if q in used_qubits:
                        continue 
                    
                    for j in graph_r.outputs():
                        if graph_r.qubit(j) == q:
                            if not graph_r.connected(i, j):
                                graph_r.add_edge((i, j))
                            break
                
                for v in graph_r.outputs():
                    if v not in graph_r.vertex_set():
                        graph_r.add_vertex_indexed(v)
                        graph_r.set_type(v, VertexType.BOUNDARY)
                        graph_r.set_row(v, graph_r.depth())
                        graph_r.set_qubit(v, self.total_graph.qubit(v))
                                        
                clean_orphan_boundaries(graph_r, temp_inputs, temp_outputs)
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
                in1,in2,out1 = create_rough_merge(graph_r, self.total_graph, q_in1, q_in2, q_out1, measurement_label)
                measurement_counter += 1
                
                temp_inputs.append(in1)
                temp_inputs.append(in2)
                temp_outputs.append(out1)
                
                non_matching = self.get_non_matching_qubit_keys(graph_r.qubits(), True)
                
                used_qubits = {q_in1, q_in2}  

                for i in non_matching:
                    q = self.total_graph.qubit(i)
                    if q in used_qubits:
                        continue  
                    
                    if i not in graph_r.vertex_set():                        
                        graph_r.add_vertex_indexed(i)
                        graph_r.set_type(i, VertexType.BOUNDARY)
                        graph_r.set_row(i, 0)
                        graph_r.set_qubit(i, q)
                        temp_inputs.append(i)

                
                uncovered_output_indices = self.get_uncovered_output_indices(graph_r)

                for i in uncovered_output_indices:
                    q = self.total_graph.qubit(i)
                    if q in used_qubits:
                        continue  
                    
                    if i not in graph_r.vertex_set():
                        graph_r.add_vertex_indexed(i)
                        graph_r.set_type(i, VertexType.BOUNDARY)
                        graph_r.set_row(i, graph_r.depth())
                        graph_r.set_qubit(i, q)
                        temp_outputs.append(i)


                temp_inputs.sort()
                temp_outputs.sort()
                
                graph_r.set_inputs((temp_inputs))
                graph_r.set_outputs((temp_outputs))
                
                for i in graph_r.inputs():
                    q = graph_r.qubit(i)
                
                    if q in used_qubits:
                        continue 
                    
                    for j in graph_r.outputs():
                        if graph_r.qubit(j) == q:
                            if not graph_r.connected(i, j):
                                graph_r.add_edge((i, j))
                            break
                
                for v in graph_r.outputs():
                    if v not in graph_r.vertex_set():
                        graph_r.add_vertex_indexed(v)
                        graph_r.set_type(v, VertexType.BOUNDARY)
                        graph_r.set_row(v, graph_r.depth())
                        graph_r.set_qubit(v, self.total_graph.qubit(v))
                                    
                clean_orphan_boundaries(graph_r, temp_inputs, temp_outputs)
                self.total_graph.compose(graph_r)
                
            else:
                raise ValueError("invalid operation")
                
    def list_operations(self):
        return self.operations

    def __repr__(self):
        return f"ParseInstructions(num_qubits={self.num_qubits}, operations={self.operations})"

# TODO: Implement this function
def initialize_ancillas(graph: Graph, qubits: list, basis: str):
    return
    
    
def max_output_index(graph: Graph):
    """Returns the maximum output index of the graph."""
    sorted_outputs = sorted(graph.outputs())

    return sorted_outputs[-1] if sorted_outputs else None

def create_smooth_split(graph: Graph, full_graph: Graph, inp: int, out1: int, out2: int, meas_outcome):
    """
    Create a smooth split (Z-spider) on `inp`, outputting to `out1` and `out2`,
    with an X correction on one leg (for Pauli frame tracking).

    Args:
        graph: PyZX graph to which new vertices are added.
        full_graph: The full graph (for boundary lookup, row alignment).
        inp: Qubit index for the input wire
        out1: Qubit index for one output (typically "upper")
        out2: Qubit index for other output (typically "lower")
        meas_outcome: A symbolic label (e.g., sympy.Symbol("m1")) for the measurement outcome

    Returns:
        Tuple of (out1_vertex, out2_vertex, input_vertex)
    """
    inp_vert = get_boundary_vertex_by_qubit(full_graph, full_graph.inputs(), inp)
    q_in = full_graph.qubit(inp_vert)
    row_base = full_graph.row(inp_vert)

    max_out = max_output_index(full_graph)

    vetices_full_graph = max_out + 1

    # input boundary vertex with same index as in full_graph
    graph.add_vertex_indexed(inp_vert)
    graph.set_type(inp_vert, VertexType.BOUNDARY)
    graph.set_qubit(inp_vert, q_in)
    graph.set_row(inp_vert, row_base)

    # Add Z-spider 
    graph.add_vertex_indexed(vetices_full_graph)
    graph.set_type(vetices_full_graph, VertexType.Z)
    graph.set_phase(vetices_full_graph, 0)
    graph.set_row(vetices_full_graph, row_base + 1)
    graph.set_qubit(vetices_full_graph, inp)
    
    split = vetices_full_graph
    
    # Add X correction spider 
    graph.add_vertex_indexed(vetices_full_graph+1)
    graph.set_type(vetices_full_graph+1, VertexType.X)
    graph.set_phase(vetices_full_graph+1, meas_outcome)
    graph.set_row(vetices_full_graph+1, row_base + 2)
    graph.set_qubit(vetices_full_graph+1, out1)

    correction = vetices_full_graph+1

    # output vertices
    graph.add_vertex_indexed(vetices_full_graph+2)
    graph.set_type(vetices_full_graph+2, VertexType.BOUNDARY)
    graph.set_qubit(vetices_full_graph+2, out1)
    graph.set_row(vetices_full_graph+2, row_base + 3)
    
    upper = vetices_full_graph+2
    
    graph.add_vertex_indexed(vetices_full_graph+3)
    graph.set_type(vetices_full_graph+3, VertexType.BOUNDARY)
    graph.set_qubit(vetices_full_graph+3, out2)
    graph.set_row(vetices_full_graph+3, row_base + 3)
    
    lower = vetices_full_graph+3

    graph.add_edges([
        (inp_vert, split),
        (split, correction),
        (correction, upper),
        (split, lower)
    ])

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
    correction = graph.add_vertex(VertexType.Z, phase=meas_outcome, qubit=out1, row=2)
    upper = graph.add_vertex(VertexType.BOUNDARY, qubit=out1, row=3)
    lower = graph.add_vertex(VertexType.BOUNDARY, qubit=out2, row=3)

    graph.add_edges([(inp_vert, split), (split, correction),
                     (correction, upper), (split, lower)])

    graph.set_inputs([inp_vert])
    graph.set_outputs([upper, lower])
    return upper, lower, inp_vert

def get_or_create_boundary(graph, vertex_id, qubit, row):
    if vertex_id in graph.vertex_set():
        return vertex_id
    graph.add_vertex_indexed(vertex_id)
    graph.set_type(vertex_id, VertexType.BOUNDARY)
    graph.set_qubit(vertex_id, qubit)
    graph.set_row(vertex_id, row)
    return vertex_id

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

    inp1 = get_boundary_vertex_by_qubit(full_graph, full_graph.inputs(), in1)
    inp2 = get_boundary_vertex_by_qubit(full_graph, full_graph.inputs(), in2)

    q_in_1 = full_graph.qubit(inp1)
    q_in_2 = full_graph.qubit(inp2)

    ## Input 1
    inp1 = get_or_create_boundary(graph, inp1, q_in_1, 0)
    graph.set_type(inp1, VertexType.BOUNDARY)
    graph.set_row(inp1, full_graph.row(q_in_1))
    graph.set_qubit(inp1, q_in_1)
    
    ## Output
    n_out = graph.add_vertex(VertexType.BOUNDARY, qubit=out, row=full_graph.row(q_in_1)+3)

    ## Input 2
    inp2 = get_or_create_boundary(graph, inp2, q_in_2, 0)
    graph.set_type(inp2, VertexType.BOUNDARY)
    graph.set_row(inp2, full_graph.row(q_in_1))
    graph.set_qubit(inp2, q_in_2)
    
    last_vertex_id = max(full_graph.vertex_set())+1
    
    
    ## Correction  correction = graph.add_vertex(VertexType.X, phase=meas_outcome*pi, qubit=in1, row=inp_row+1)
    graph.add_vertex_indexed(last_vertex_id)
    graph.set_type(last_vertex_id, VertexType.X)
    graph.set_phase(last_vertex_id, meas_outcome)
    graph.set_row(last_vertex_id, full_graph.row(q_in_1)+1)
    graph.set_qubit(last_vertex_id, in1)
    
    ## Merge merge = graph.add_vertex(VertexType.Z, phase=0, qubit=in1, row=inp_row+2)
    graph.add_vertex_indexed(last_vertex_id+1)
    graph.set_type(last_vertex_id+1, VertexType.Z)
    graph.set_phase(last_vertex_id+1, 0)
    graph.set_row(last_vertex_id+1, full_graph.row(q_in_1)+2)
    graph.set_qubit(last_vertex_id+1, in1)
    

    graph.add_edges([(inp1, last_vertex_id), (last_vertex_id, last_vertex_id+1),
                     (inp2, last_vertex_id+1), (last_vertex_id+1, n_out)])

    graph.set_inputs([in1, in2])
    graph.set_outputs([n_out])

    # If we want to track the correction in a PauliFrame:
    if pauli_frame is not None:
        pauli_frame.update(in1, "X")  # record X byproduct on qubit in1

    return inp1, inp2, n_out


def create_rough_merge(graph: Graph, full_graph: Graph, in1: int, in2: int, out: int, meas_outcome, pauli_frame: PauliFrame=None):
    """Create a smooth merge (Z-spider) with X correction on upper input leg.
    
    Args:
        graph: ZX graph to build upon
        in1: Qubit index for upper input boundary
        in2: Qubit index for lower input boundary
        out: Qubit index for output boundary
        meas_outcome: Symbolic measurement outcome
        pauli_frame: Optional PauliFrame object to record the X byproduct
    """

    inp1 = get_boundary_vertex_by_qubit(full_graph, full_graph.inputs(), in1)
    inp2 = get_boundary_vertex_by_qubit(full_graph, full_graph.inputs(), in2)

    q_in_1 = full_graph.qubit(inp1)
    q_in_2 = full_graph.qubit(inp2)

    ## Input 1
    inp1 = get_or_create_boundary(graph, inp1, q_in_1, 0)
    graph.set_type(inp1, VertexType.BOUNDARY)
    graph.set_row(inp1, full_graph.row(q_in_1))
    graph.set_qubit(inp1, q_in_1)
    
    ## Output
    n_out = graph.add_vertex(VertexType.BOUNDARY, qubit=out, row=full_graph.row(q_in_1)+3)

    ## Input 2
    inp2 = get_or_create_boundary(graph, inp2, q_in_2, 0)
    graph.set_type(inp2, VertexType.BOUNDARY)
    graph.set_row(inp2, full_graph.row(q_in_1))
    graph.set_qubit(inp2, q_in_2)
    
    last_vertex_id = max(full_graph.vertex_set())+1
    
    
    ## Correction  correction = graph.add_vertex(VertexType.X, phase=meas_outcome*pi, qubit=in1, row=inp_row+1)
    graph.add_vertex_indexed(last_vertex_id)
    graph.set_type(last_vertex_id, VertexType.Z)
    graph.set_phase(last_vertex_id, meas_outcome)
    graph.set_row(last_vertex_id, full_graph.row(q_in_1)+1)
    graph.set_qubit(last_vertex_id, in1)
    
    ## Merge merge = graph.add_vertex(VertexType.Z, phase=0, qubit=in1, row=inp_row+2)
    graph.add_vertex_indexed(last_vertex_id+1)
    graph.set_type(last_vertex_id+1, VertexType.X)
    graph.set_phase(last_vertex_id+1, 0)
    graph.set_row(last_vertex_id+1, full_graph.row(q_in_1)+2)
    graph.set_qubit(last_vertex_id+1, in1)
    

    graph.add_edges([(inp1, last_vertex_id), (last_vertex_id, last_vertex_id+1),
                     (inp2, last_vertex_id+1), (last_vertex_id+1, n_out)])

    graph.set_inputs([in1, in2])
    graph.set_outputs([n_out])

    # If we want to track the correction in a PauliFrame:
    if pauli_frame is not None:
        pauli_frame.update(in1, "X")  # record X byproduct on qubit in1

    return inp1, inp2, n_out



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