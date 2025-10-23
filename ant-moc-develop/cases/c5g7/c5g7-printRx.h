#ifdef ENABLE_MPI_
  if (mpigod::isMPIRoot()) {
#endif
    log_printf(TITLE, "Reaction Rates");
    Mesh mesh(solver);
    mesh.createLattice(51, 51, 45);
    for (int i = -1; i < 8; ++i) {
      mesh.printMeshDataToXML("rx", "Fission reaction rates...", TallyType::Fission_RX, false, i);
      mesh.printMeshDataToXML("rx", "Scalar fluxes...", TallyType::Flux_RX, false, i);
    }

    log_printf(TITLE, "Volume-averaged Reaction Rates");
    for (int i = -1; i < 8; ++i) {
      mesh.printMeshDataToXML("rx", "Fission reaction rates(volume-averaged)...", TallyType::Fission_RX, true, i);
      mesh.printMeshDataToXML("rx", "Scalar fluxes(volume-averaged)...", TallyType::Flux_RX, true, i);
    }


#ifdef ENABLE_MPI_
  }
#endif
