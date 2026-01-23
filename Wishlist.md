# QuLTRA – Wishlist

This document collects the main features, improvements, and development goals planned for **QuLTRA**.

---

## 1. More Testing & Validation

- Thoroughly test **QuLTRA as a standalone tool** using different circuit configurations.
- Include tests with:
  - Different circuit topologies
  - A large number of circuit elements
- Verify numerical stability, scalability, and performance as circuit complexity increases.

---

## 2. Integration with Qiskit Metal

Investigate and validate a workflow to interface **QuLTRA** with **Qiskit Metal**, in particular:

- Create circuit layouts using **Qiskit Metal**
- Extract the **capacitance matrix** using electromagnetic solvers:
  - Q3D
  - Similar EM tools
- Use the extracted capacitance matrix as input to build the corresponding circuit in **QuLTRA**
- Complete the circuit by adding missing elements, such as:
  - Josephson junctions
  - Lumped inductors
  - CPWs
  - CPW couplers
- Verify that the full workflow functions correctly and produces consistent results.

---

## 3. Code Refactoring & Organization

- Reorganize the codebase to improve readability and maintainability
- Assign **one Python file per class**
- Clarify module structure and dependencies
- Improve overall code cleanliness and documentation

---

## 4. Graphical User Interface (GUI)

- Develop a **simple GUI** that allows users to:
  - Build and edit circuits
  - Define circuit elements and connections
- The initial goal is usability rather than full feature coverage

---

This wishlist will evolve as the project grows and new requirements emerge.
