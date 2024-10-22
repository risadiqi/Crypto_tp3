# TP3 RSA

Dans ce TP, nous avons implémenté l'algorithme de chiffrement RSA en utilisant la bibliothèque GMP. RSA est l'un des algorithmes de cryptographie asymétrique les plus utilisés au monde. Il repose sur la difficulté de factoriser de grands nombres premier. RSA se base sur la génération d'une paire de clés : une clé publique pour le chiffrement et une clé privée pour le déchiffrement. Pour générer les clés, on choisit deux nombres premiers p et q, puis on calcule n = p * q et le qotient d'Euler ϕ(n) = (p-1) * (q-1). Un exposant de chiffrement e est choisi tel qu'il soit premier avec ϕ(n), puis l'inverse modulaire d est calculé pour le déchiffrement. Le chiffrement s'effectue en élevant le message à la puissance e modulo n, et le déchiffrement consiste à élever le message chiffré à la puissance d modulo n. La sécurité de RSA repose sur la difficulté de factoriser n en p et q.

## Membres de l'équipe.

 + Nouhaila Jabbar
 + Rim Sadiqi

## Détails des fichiers

* Squelette de Code/ : Ce répertoire contient le code source principal :

     + main.cpp : Ce fichier implémente l'algorithme RSA en passant par toutes les étapes de chiffrement et de déchiffrement.
 
     + main_test.cpp : Ce fichier n'a pas été utilisé.
       
 



