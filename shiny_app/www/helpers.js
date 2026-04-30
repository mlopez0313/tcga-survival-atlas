// helpers.js — small tweaks for the TCGA Survival Atlas Shiny app.
// Currently a placeholder for future client-side niceties (keyboard shortcuts,
// tooltip wiring, etc.). Intentionally minimal so it never blocks app load.

document.addEventListener("DOMContentLoaded", function () {
  // Smooth-scroll inside any nav links pointing at internal anchors.
  document.querySelectorAll('a[href^="#"]').forEach(function (a) {
    a.addEventListener("click", function (e) {
      const target = document.querySelector(a.getAttribute("href"));
      if (target) {
        e.preventDefault();
        target.scrollIntoView({ behavior: "smooth", block: "start" });
      }
    });
  });
});
