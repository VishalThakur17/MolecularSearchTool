(function () {
    function parseJsonDataset(element, key, fallback) {
        try {
            return JSON.parse(element.dataset[key] || "null") || fallback;
        } catch (error) {
            console.warn("Could not parse dataset:", key, error);
            return fallback;
        }
    }

    function rawStructureUrl(structure) {
        if (!structure) return "";

        const sourceType = (structure.source_type || "").toLowerCase();
        const sourceId = (structure.source_id || "").trim();

        if (!sourceId) return "";

        if (sourceType === "pdb") {
            return `https://files.rcsb.org/download/${sourceId.toUpperCase()}.pdb`;
        }

        if (sourceType === "afdb") {
            return `https://alphafold.ebi.ac.uk/files/AF-${sourceId}-F1-model_v4.pdb`;
        }

        return "";
    }

    function normalizeColor(hexColor, fallback) {
        const value = (hexColor || fallback || "#dc2626").replace("#", "");
        return parseInt(value, 16);
    }

    function escapeHtml(value) {
        return String(value || "")
            .replace(/&/g, "&amp;")
            .replace(/</g, "&lt;")
            .replace(/>/g, "&gt;")
            .replace(/"/g, "&quot;")
            .replace(/'/g, "&#039;");
    }

    function annotationLabel(annotation) {
        const binder = annotation.binder_name || "Binder";
        const region = `${annotation.region_start}-${annotation.region_end}`;
        const evidence = annotation.evidence_type || "Binding region";
        return `${binder}: ${region} (${evidence})`;
    }

    function buildAnnotationList(container, annotations) {
        const list = container.querySelector("[data-annotation-list]");
        if (!list) return;

        if (!annotations || annotations.length === 0) {
            list.innerHTML = `<p class="empty">No binding-site annotations are available for this 3D view.</p>`;
            return;
        }

        list.innerHTML = annotations.map((ann) => {
            const color = ann.color_hex || "#dc2626";
            const inferredText = ann.inferred ? "Inferred" : "Curated";
            const sourceNote = ann.source_note || "No source note available.";

            return `
                <div class="three-d-annotation-row">
                    <span class="three-d-dot" style="background:${escapeHtml(color)};"></span>
                    <div>
                        <strong>${escapeHtml(ann.binder_name || "Unknown binder")}</strong>
                        <p>
                            Residues ${escapeHtml(ann.region_start)}-${escapeHtml(ann.region_end)}
                            · ${escapeHtml(ann.region_label || "Binding region")}
                            · ${escapeHtml(inferredText)}
                        </p>
                        <p class="three-d-source">${escapeHtml(sourceNote)}</p>
                    </div>
                </div>
            `;
        }).join("");
    }

    function showStatus(container, message, type) {
        const status = container.querySelector("[data-structure-status]");
        if (!status) return;

        status.textContent = message || "";
        status.className = `three-d-status ${type || ""}`;
    }

    async function loadStructure(container, viewer, structure, annotations) {
        const url = rawStructureUrl(structure);

        viewer.clear();

        if (!url) {
            showStatus(container, "No usable PDB or AlphaFold structure URL is available.", "error");
            viewer.render();
            return;
        }

        showStatus(container, "Loading 3D structure...", "loading");

        try {
            const response = await fetch(url);

            if (!response.ok) {
                throw new Error(`Structure download failed: ${response.status}`);
            }

            const pdbText = await response.text();

            viewer.addModel(pdbText, "pdb");

            viewer.setStyle({}, {
                cartoon: {
                    color: "lightgray",
                    opacity: 0.85
                }
            });

            annotations.forEach((ann, index) => {
                const start = Number(ann.region_start);
                const end = Number(ann.region_end);

                if (!Number.isFinite(start) || !Number.isFinite(end)) return;

                const color = normalizeColor(ann.color_hex, "#dc2626");

                viewer.setStyle(
                    { resi: `${start}-${end}` },
                    {
                        cartoon: {
                            color: color,
                            opacity: 1.0
                        },
                        stick: {
                            color: color,
                            radius: 0.22
                        }
                    }
                );

                viewer.addLabel(
                    annotationLabel(ann),
                    {
                        position: { x: 0, y: index * 2.2, z: 0 },
                        backgroundColor: "white",
                        fontColor: "black",
                        borderColor: "gray",
                        borderThickness: 1,
                        fontSize: 11,
                        inFront: false
                    }
                );
            });

            viewer.zoomTo();
            viewer.render();

            showStatus(
                container,
                `Loaded ${structure.label || structure.source_id}. Highlighted ${annotations.length} binding region(s).`,
                "success"
            );

            buildAnnotationList(container, annotations);

        } catch (error) {
            console.error(error);
            showStatus(
                container,
                "Unable to load this 3D structure. The 2D binding map still remains available.",
                "error"
            );
            viewer.render();
        }
    }

    function init3DHighlightViewer(container) {
        const viewerElement = container.querySelector("[data-structure-viewer]");
        const select = container.querySelector("[data-structure-select]");

        if (!viewerElement || !window.$3Dmol) return;

        const structures = parseJsonDataset(container, "structures", []);
        const annotations = parseJsonDataset(container, "annotations", []);

        if (!structures || structures.length === 0) {
            showStatus(container, "No 3D structure is available for this page.", "error");
            return;
        }

        viewerElement.style.position = "relative";
        viewerElement.style.overflow = "hidden";
        viewerElement.innerHTML = "";

        const viewer = window.$3Dmol.createViewer(viewerElement, {
            backgroundColor: "white"
        });

        if (select) {
            select.innerHTML = structures.map((item, index) => {
                const label = item.label || item.source_id || `Structure ${index + 1}`;
                return `<option value="${index}">${escapeHtml(label)}</option>`;
            }).join("");

            select.addEventListener("change", function () {
                const selectedStructure = structures[Number(select.value)] || structures[0];
                loadStructure(container, viewer, selectedStructure, annotations);
            });
        }

        loadStructure(container, viewer, structures[0], annotations);
    }

    document.addEventListener("DOMContentLoaded", function () {
        document.querySelectorAll("[data-binding-3d-viewer]").forEach(init3DHighlightViewer);
    });
})();
