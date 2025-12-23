# Pour la compilation et l'execution 🎀



## Gauss Seidel

Pour compiler et exécuter Gauss-Seidel en statique faire:

```bash
cd GS-static
gcc -Wall main.c gs_static.c utilities.c -o exec -lm
./exec
```
Pour compiler et exécuter Gauss-Seidel en dynamique faire:

```bash
cd GS-dynamic
gcc -Wall main.c gs_dynamic.c utilities.c blas1_static.c -o exec -lm
./exec
```

## SOR

Pour compiler et exécuter SOR en statique faire:
```bash
cd SOR-static
gcc -Wall main.c sor_static.c utilities.c blas1_static.c -o exec -lm
./exec
```






