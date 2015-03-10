        GlobalBasis global_basis;

        GlobalBasis::LocalView local_view(global_basis);

        local_view.globalBasis();
        local_view.bind(e);
        size_type N_TH = local_view.size();
        size_type N_TH_max = local_view.maxSize();
        LocalView::LeafNode pressure = local_view.tree().child(...);
        size_type N_p = pressure.size();
        pressure.localIndex(i); // [0..N_p) -> [0..N_TH)
        auto& fe = pressure.finiteElement();

        // Variante 5
        auto index_set = global_basis.indexSet();
        index_set.size(prefix); // TBD
        auto local_index_set = index_set.localIndexSet();
        local_index_set.bind(local_view); // berechnet globale Indizes
        local_index_set.index(i);
        local_index_set.size(); // == localView().size()
        local_index_set.localView();
        local_index_set.unbind();
