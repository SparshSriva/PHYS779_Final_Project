# ZXLS.py

from pyzx import Graph, VertexType, simplify
from sympy import pi

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