# Pytania i uwagi do metaanalizy

## 1. PROBASE — problem klasyfikacji

### Główna obserwacja

PROBASE nie pasuje do żadnej z dwóch zdefiniowanych rodzin:

| Kategoria | Definicja / rzeczywisty kontrast |
|---|---|
| **Rodzina A** | Screening invitation vs no screening / usual care |
| **Rodzina B** | Diagnostic / biopsy strategy (both arms screened) |
| **PROBASE** | I: PSA od razu w wieku 45 lat; C: PSA odroczone do 50. r.ż. |

- **Nie pasuje do rodziny B**, ponieważ w momencie analizy (runda 1) ramię C nie miało PSA. Narusza to definicję „both arms screened”.
- **Nie pasuje do rodziny A**, ponieważ C nie oznacza „no screening”. Jest to ta sama interwencja przesunięta w czasie. Oba ramiona zaproszono na screening; randomizacja dotyczyła momentu rozpoczęcia, a nie samego faktu screeningu.

To trzeci, osobny kontrast: **timing of screening initiation**, którego nie uwzględnia obecny schemat.

### Co faktycznie działo się w ramieniu C?

Nie było to ramię bez diagnostyki:

- `bx = 37`,
- `pca = 2`,
- podgrupa 6537 mężczyzn otrzymała DRE, co doprowadziło do 37 biopsji i wykrycia 2 raków.

PROBASE-C było więc ramieniem z DRE w podgrupie, a nie ramieniem całkowicie nieprzebadanym.

### Odstępstwo od protokołu

Protokół (§113) wymienia grupy strategii, a pierwszą z nich jest samodzielne DRE:

> “including digital rectal examination (DRE), PSA-based screening, PSA combined with DRE, PSA followed by MRI…”

W analizie nie ma osobnej grupy DRE. Jedyne ramię, które mogłoby do niej trafić — PROBASE-C — otrzymało etykietę PSA. Ten sam problem występuje więc na dwóch poziomach: przypisania do strategii oraz przypisania do rodziny.

### Dlaczego `RR = 7.487` nie jest porównywalne z pozostałymi?

Działają tu dwa niezależne mechanizmy:

1. **Kontrast dotyczy czasu, a nie samego screeningu.**  
   Porównanie brzmi „przebadaj teraz” vs „przebadaj za 5 lat”, ale outcome odczytano przed drugim screeningiem. Wynik jest zatem w dużej mierze zdefiniowany przez konstrukcję badania: jedna grupa przeszła diagnostykę, a druga jeszcze nie.

2. **Wiek 45 lat i jedna runda oznaczają niemal zerową zachorowalność tła.**  
   W ERSPC, CAP i ProScreen ramię kontrolne przez 3–10 lat akumuluje raki wykryte klinicznie, więc RR mierzy nadwyżkę wykrycia netto — sensowną miarę overdiagnosis. U 45-latków obserwowanych przez jedną rundę rozpoznań klinicznych jest praktycznie zero. Dwa raki pochodziły z biopsji po DRE, a nie z diagnostyki objawowej. Mianownik jest sztucznie bliski zeru, dlatego `RR = 7.487` odzwierciedla punkt odczytu i wiek kohorty, a nie to samo zjawisko co w pozostałych badaniach.

### Wpływ liczbowy

PROBASE ma tylko 3,4% wagi (17 zdarzeń), więc wpływ na dwa pierwsze outcomes jest mały, ale na trzeci — decydujący:

| Outcome w rodzinie A | Z PROBASE | Bez PROBASE |
|---|---:|---:|
| GG1 detection | 1.756 (1.085–2.842), `p = 0.031` | 1.668 (1.118–2.488), `p = 0.027` |
| csPCa detection | 1.156 (0.756–1.768) | 1.138 (0.755–1.716) |
| Grade mix (GG1/PCa) | 1.079 (0.639–1.823) | 1.210 (0.998–1.468) |
| Test różnicy A vs B (grade mix) | `p = 0.019` | `p = 0.00025` |

Usunięcie PROBASE nieco poprawia dwa pierwsze wyniki. W przypadku grade mix — głównego wyniku pracy — PROBASE istotnie zmienia rezultat:

- ramię kontrolne ma `2/2 = 100% GG1`,
- własne `RR = 0.320`,
- pooled estimate rodziny A spada z 1.21 do 1.08 i traci istotność.

Wartość 100% przy `n = 2` nie daje wiarygodnego oszacowania.

### Wniosek

Obecny raport traktuje wykluczenie PROBASE jako sensitivity analysis (§3.2), a wersję z PROBASE jako analizę główną. Jeśli PROBASE nie spełnia definicji rodziny A, hierarchia powinna być odwrotna:

- **analiza główna:** bez PROBASE,
- **sensitivity analysis:** z PROBASE.

Główny wynik wynosiłby wtedy:

> **1.210 (0.998–1.468), `p = 0.00025` dla różnicy między rodzinami.**

To zmiana wymagająca uzgodnienia z Miszczykiem i Rajwą, ponieważ dotyczy definicji analizy głównej względem zarejestrowanego protokołu.

Możliwe rozwiązania:

1. Wykluczyć PROBASE z syntezy porównawczej i opisać narracyjnie.
2. Utworzyć osobną rodzinę „timing” z `k = 1` — praktycznie równoważne rozwiązaniu 1.
3. Pozostawić PROBASE w rodzinie A, ale jasno opisać w ograniczeniach, że kontrast i punkt odczytu różnią się od pozostałych czterech badań.

**Do rozważenia:** przeliczenie wariantu, w którym PROBASE jest wykluczone z całej Part 1 i Part 2, a nie tylko z RR w §3.2.

---

## 2. Dlaczego `p_C` w ProScreen jest o rząd wielkości mniejsze?

> Dlaczego `p_C` w ProScreen jest o rząd wielkości mniejsze niż w CAP, ERSPC i PLCO, skoro w każdym z tych badań ramię kontrolne oznacza brak screeningu?

### Krótka odpowiedź

„Brak screeningu” nie oznacza „braku diagnostyki”. Ponadto wartości odnoszą się do różnych okresów obserwacji. Łączny, około 33-krotny rozrzut można rozłożyć na trzy mnożniki po około 3×.

> Liczby poniżej przeliczono ręcznie na podstawie wartości z `results.txt` i `meta-analysis.csv`.

### Dekompozycja

| Ramię kontrolne | Kraj / region | Lata obserwacji | GG1 surowo | GG1/rok | any PCa/rok | Udział GG1 |
|---|---|---:|---:|---:|---:|---:|
| PLCO-C | USA | 10.0 | 4.675% | 0.468 | 0.775 | 60.3% |
| ERSPC-C | Europa | 9.0 | 2.641% | 0.293 | 0.536 | 54.8% |
| CAP-C | UK | 10.0 | 1.112% | 0.111 | 0.358 | 31.1% |
| ProScreen-C | Finlandia | 3.2 | 0.143% | 0.045 | 0.244 | 18.3% |

```text
PLCO / ProScreen = 32.8×
                  ≈ 3.1×  czas obserwacji (10 lat vs 3.2 roku)
                  × 3.2×  zapadalność na jakikolwiek PCa na rok
                  × 3.3×  odsetek wykrytych raków sklasyfikowanych jako GG1
```

#### Mnożnik 1 — czas obserwacji

To kumulacyjna incydencja, nie prewalencja. W nieprzebadanej populacji liczba ujawnionych raków rośnie wraz z czasem obserwacji. Różnica 10 lat vs 3,2 roku daje około 3×.

#### Mnożnik 2 — oportunistyczny PSA

W ramieniu bez zaplanowanego screeningu bezobjawowy GG1 zostanie znaleziony głównie wtedy, gdy ktoś wykona PSA poza badaniem. `p_C` mierzy więc w dużej mierze intensywność oportunistycznego („dzikiego”) screeningu w danym kraju i epoce, a nie wyłącznie biologię:

- **PLCO (USA, 1993–2001):** ekstrakcja podaje 40–52% PSA w ramieniu kontrolnym — najwyższa wartość;
- **CAP (UK, 2001–2009):** NHS praktycznie nie wykonywał PSA — wartość roczna około 2× niższa niż w PLCO;
- **ProScreen (Finlandia, 2018+):** krótka obserwacja i współczesna praktyka diagnostyczna.

#### Mnożnik 3 — klasyfikacja GG1

Udział GG1 wśród wykrytych raków spada:

> **60% → 55% → 31% → 18%**

Dzieje się to w ramionach bez badanej interwencji. Nie może to więc wynikać wyłącznie ze strategii screeningu.

### Dlaczego nie psuje to RR?

To uzasadnienie, dlaczego Part 1 pozostaje solidna:

| Trial | `p_I` | `p_C` | RR |
|---|---:|---:|---:|
| ERSPC | 5.934% | 2.641% | 2.247 |
| PLCO | 5.918% | 4.675% | 1.266 |
| CAP | 1.723% | 1.112% | 1.550 |
| ProScreen | 0.257% | 0.143% | 1.798 |

- rozrzut `p_C`: **32.8×**,
- rozrzut `p_I`: **23.1×**,
- rozrzut RR: **1.8×**.

Czas obserwacji, wiek, oportunistyczny PSA i konwencja gradingu działają podobnie na oba ramiona tego samego badania, dlatego w dużej mierze skracają się w ilorazie. RR przy surowych stopach rozrzuconych 33-krotnie pozostaje w zakresie 1.27–2.25.

To właśnie oznacza „ochronę przez randomizację”: nie tylko przed confoundingiem osobniczym, ale również przed czynnikami wspólnymi dla obu ramion.

Dlatego `τ² = 0.084` (SD 0.29 na skali log) jest umiarkowane, mimo że surowe dane różnią się o rząd wielkości.

### Co ten problem faktycznie podważa?

#### 1. Grupę „No screening / usual care” w Part 2

Pooled estimate wynosi 1.20% (0.33–4.30), przy `I² = 99.5%`. Uśrednia cztery ramiona różniące się 33-krotnie z przyczyn niezwiązanych bezpośrednio ze strategią. Nie jest to jedna interpretowalna estymata, lecz średnia z czterech różnych kontekstów.

#### 2. Porównywalność gradingu między epokami

Migracja klasyfikacji Gleasona zagraża głównemu wynikowi pracy, a `REPORT.md` obecnie jej nie omawia.

Rewizje ISUP z 2005 i 2014 przesunęły wzory cribriform i część wzoru 3 do wzoru 4. W efekcie guz klasyfikowany w latach 90. jako Gleason 6 może być dziś sklasyfikowany jako 3+4, czyli GG2. Dodatkowo kategorie Gleason 2–4 przestały istnieć, podczas gdy ekstrakcja PLCO definiuje GG1 jako „Gleason 2–4 + Gleason 5–6”.

Protokół (§6) mówi:

> “When ISUP grading was not explicitly reported, Gleason score ≤6 will be considered equivalent to ISUP grade group 1.”

Ta równoważność jest definicyjnie prawdziwa, ale nie zapewnia porównywalności między epokami. Ten sam guz może otrzymać inny wynik zależnie od okresu, w którym oceniał go patolog.

Ramiona kontrolne działają tu jak negative control:

```text
Gradient „ery”       (ramiona bez interwencji): 60.3% → 18.3% = 3.3×
Gradient „strategii” (teza pracy, §2.3):        65.7% → 16.4% = 4.0×
```

Rzekomy efekt strategii jest tego samego rzędu co efekt kalendarza mierzony w ramionach bez interwencji. Era i strategia są w tych danych silnie skorelowane (`PSA + DRE = 2009`, `biomarker → MRI = 2021–2024`), więc nie da się ich rozdzielić w analizie przekrojowej.

Nie oznacza to, że teza jest fałszywa. Oznacza natomiast, że §2.3 i §5 pkt 2 mierzą łącznie:

1. **Realny efekt triage’u** — mniej niegroźnych raków trafia do biopsji.
2. **Migrację gradingu** — ten sam guz jest punktowany wyżej, co stanowi artefakt pomiarowy.
3. **Biopsję targeted vs systematic** — biopsja celowana trafia w zmianę indeksową, a systematyczna wykrywa przypadkowe mikroogniska GG1; ten element jest uzasadnioną częścią strategii.

### Rekomendowane działania

Wersja tezy oparta na RR wewnątrz badania jest odporna na migrację gradingu, ponieważ oba ramiona były oceniane przez patologów w tym samym okresie. Dotyczy to rodziny B:

> **Grade mix `RR = 0.567 (0.302–1.064)`**

To ten wynik powinien być rozważony jako główny, zamiast tabeli z §2.3.

1. **Dopisać migrację gradingu jako ograniczenie.**  
   Podać gradient ramion kontrolnych (60.3% / 54.8% / 31.1% / 18.3%) jako jego ilościową miarę.

2. **Przeprowadzić analizę wrażliwości ograniczoną do jednej epoki gradingu (post-2014).**  
   Uwzględnić: ProScreen, STHLM3-MRI, MVP, GÖTEBORG-2, PROSA i PROBASE.

   Po takim ograniczeniu:

   - PSA + DRE (65.7%) znika całkowicie — było to jedno ramię, PLCO-I z 2009 roku;
   - grupa PSA obejmuje tylko małe ramiona: PROBASE-I 31.2%, STHLM3-C 40.8%, MVP-C 50.0% oraz problematyczne PROBASE-C 2/2;
   - biomarker → MRI (15.0% / 17.8%) pozostaje najniższy;
   - dla MRI-first (26.7% / 34.6%) przewaga nad PSA → MRI (36.4–53.4%) staje się niewyraźna.

Po ograniczeniu wpływu ery pozostaje jeden mocny wniosek:

> **Biomarker → MRI daje około 16% GG1 vs 31–53% dla PSA i PSA → MRI w tej samej epoce.**

Wynik ten jest potwierdzony wewnątrz badania przez STHLM3-MRI (`grade mix RR = 0.369`).

### Anomalia do sprawdzenia

W ProScreen:

> ProScreen-I 17.8% vs ProScreen-C 18.3% → `grade mix RR = 0.973 (0.679–1.393)`

W tym badaniu ścieżka biomarker → MRI dała praktycznie taki sam miks stopni jak brak screeningu w tej samej epoce i kraju. Efekt rodziny B opiera się więc na STHLM3-MRI i GÖTEBORG-2, a nie na ProScreen. Warto zaznaczyć to wprost, ponieważ recenzent prawdopodobnie zwróci na to uwagę.

**Do rozważenia:** dodanie do skryptu wariantu „jedna epoka gradingu”, aby sprawdzić, jaka część gradientu pozostaje po tym ograniczeniu.

----

Bardzo dobre pytanie — to jest fundament Part 2 i źródło większości problemów, które wypisałem na początku.

Najprościej: to inne GROUP BY
Masz jedną tabelę 18 wierszy (arm_level_data.csv). Part 1 i Part 2 robią na niej dwa różne zapytania:

-- PART 1
SELECT ... FROM ramiona
JOIN po badaniu (I z C)         -- jednostka = BADANIE (para ramion)
GROUP BY family                 -- 2 grupy: A i B

-- PART 2
SELECT ... FROM ramiona         -- jednostka = RAMIĘ (bez JOIN-a!)
GROUP BY strategy               -- 6 grup
W Part 1 ramiona żyją w parach — ERSPC-I istnieje tylko po to, żeby porównać je z ERSPC-C. W Part 2 pary są rozerwane, każde ramię staje się samodzielnym wierszem, i wiersze są przetasowane według zupełnie innego kryterium.

Kryterium: „co decyduje, kto idzie na biopsję"
Autorzy nadali każdemu z 18 ramion etykietę odpowiadającą na pytanie: jaki test jest filtrem, który kwalifikuje mężczyznę do biopsji? W kodzie jest to opisane jako „arms are classified by the triage test that triggers work-up".

Sześć możliwych etykiet, uszeregowanych jak epoki historyczne:

No screening  →  PSA + DRE  →  PSA  →  PSA→MRI  →  Biomarker→MRI  →  MRI-first
   (nic)        (palec+krew)   (krew)  (krew,potem  (krew+algorytm,  (MRI od razu)
                                        MRI)         potem MRI)
To jest cała teza pracy: im dalej w prawo, tym więcej filtrów przed biopsją, tym mniej wykrytego GG1.

Jak wygląda przetasowanie
Kolumna po lewej to badania (jak w Part 1), po prawej — do której grupy trafiło każde ramię:

badanie	ramię I →	ramię C →
ERSPC	PSA	No screening
PLCO	PSA + DRE	No screening
CAP	PSA	No screening
ProScreen	Biomarker → MRI	No screening
PROBASE	PSA	PSA ⚠️
STHLM3-MRI	Biomarker → MRI	PSA
MVP	MRI-first	PSA
GÖTEBORG-2	PSA → MRI	PSA → MRI ⚠️
PROSA	MRI-first	PSA → MRI
Zauważ dwie rzeczy:

ramiona jednego badania trafiają zwykle do różnych grup (ERSPC-I do „PSA", ERSPC-C do „No screening")
ale przy PROBASE i GÖTEBORG-2 oba ramiona trafiają do tej samej grupy ⚠️ — to są moje „bliźniaki" z poprzedniej odpowiedzi
Grupa „PSA" — 6 ramion z sześciu różnych miejsc
Wyciągając z tabeli wszystko z etykietą „PSA":

ramię	skąd	GG1 na przebadanego
ERSPC-I	ramię interwencyjne ERSPC	7.2176%
CAP-I	ramię interwencyjne CAP	4.8475%
MVP-C	ramię kontrolne MVP	1.6129%
STHLM3-MRI-C	ramię kontrolne STHLM3	1.4219%
PROBASE-I	ramię interwencyjne PROBASE	0.0644%
PROBASE-C	ramię kontrolne PROBASE	0.0086%
To jest kluczowa obserwacja: trzy z nich to ramiona interwencyjne, trzy to kontrolne. W Part 1 to były przeciwne strony tej samej randomizacji. W Part 2 lądują w jednym worku, bo obie strony używały PSA jako filtra.

To samo w sobie nie jest błędem — z punktu widzenia pytania „ile GG1 wykrywa strategia PSA" nie ma znaczenia, czy dane ramię było w swoim badaniu „interwencją" czy „kontrolą". Liczy się, co faktycznie robiono z pacjentami.

Problem: etykieta „PSA" to jedyne, co je łączy
ERSPC-I	CAP-I	STHLM3-C	MVP-C	PROBASE-I	PROBASE-C
wiek	55–69	50–69	50–74	~50–75	45	45
liczba rund PSA	wiele (co 4 lata)	1	1	1	1	0 ⚠️
okres zbierania raków	9 lat	10 lat	1 runda	1 runda	1 runda	1 runda
ilu faktycznie przebadano	82%	35.5%	100%	93%	100%	— ⚠️
próg PSA	≥3	≥3	≥3	≥2.6	≥3	—
n przebadanych	59 923	67 313	5 134	248	23 301	—
Porównaj skrajności:

ERSPC-I: mężczyźni 55–69 lat, PSA co 4 lata, raki zbierane przez 9 lat → 7.22%
PROBASE-I: mężczyźni 45-letni, jedno badanie PSA, raki z jednej rundy → 0.06%
Różnica 112× nie ma nic wspólnego ze „strategią PSA". Ma wszystko wspólnego z tym, że 45-latek badany raz to nie to samo co 60-latek badany pięć razy przez dekadę. Etykieta „PSA" jest identyczna, eksperyment zupełnie inny.

Stąd 841× rozrzut w grupie, tau² = 5.90 i średnia geometryczna 0.60%, która leży poniżej czterech z sześciu ramion. Artefakt nie jest błędem matematyki — jest konsekwencją tego, że grupa nie jest grupą.

I stąd dwa konkretne błędy z mojej listy
A4 — PROBASE-C nie powinno być w tej grupie. To ramię opóźnionego screeningu: ci mężczyźni w rundzie 1 nie dostali PSA wcale (tylko podgrupa 6537 miała DRE). Etykieta mówi „strategia = PSA", a faktyczna strategia to „nie rób nic jeszcze". Sprawdziłem konsekwencje:

grupa PSA jak jest (k=6)     0.60%   tau² = 5.90
bez PROBASE-C (k=5)          1.39%   tau² = 2.84
bez PROBASE w ogóle (k=4)    3.15%   tau² = 0.51
Jedna błędna etykieta zmienia wynik pięciokrotnie i zamienia grupę beznadziejnie heterogeniczną w prawie jednorodną.

Jest tu też jawna sprzeczność w danych: PROBASE ma family = A („invitation vs no screening"), ale jego ramię kontrolne dostaje strategy = PSA. Ten sam plik mówi o tym samym ramieniu dwie różne rzeczy.

A6 — bliźniaki. PROBASE-I i PROBASE-C to ci sami pacjenci, kraj, laboratorium, patolog, ta sama runda. Model traktuje je jako dwie niezależne obserwacje. Grupa ma nominalnie k=6, efektywnie ~5.

Podsumowując
„Grupa PSA" = 6 ramion z 5 różnych badań, które autorzy uznali za stosujące tę samą strategię diagnostyczną, wyrwane ze swoich par randomizacyjnych i zsypane razem.

I to jest źródło wszystkich ograniczeń Part 2:

utrata randomizacji — porównujesz ramiona z różnych badań, więc żadna liczba tutaj nie jest kauzalna
etykieta ≠ eksperyment — dwa ramiona z tą samą etykietą mogą różnić się wiekiem, liczbą rund i długością obserwacji
klasyfikacja jest decyzją autorów, nie faktem z danych — i jak pokazuje PROBASE-C, jedna zmiana etykiety przesuwa wynik 5×
Dlatego Part 2 jest opisowy, a Part 1 kauzalny — i dlatego weryfikacja Part 2 polega głównie na sprawdzeniu etykiet, nie matematyki.

---

Zajrzałem do protocol.md. Odpowiedź jest jednoznaczna i niewygodna.

Co protokół faktycznie mówi
§6 Endpoint Definitions, primary endpoint — dosłownie:

„Detection rate and proportion of ISUP grade group 1 prostate cancer among screened men or among detected prostate cancers."

§9 Data Synthesis — jeszcze bardziej wprost:

„The proportion of ISUP grade group 1 and ISUP grade group ≥2 cancers will be presented in relation to the total number of men undergoing screening within each study arm."

i lista wyliczeniowa:

„– Detection rates per number of screened men, where denominators are available.
– The proportion of ISUP grade group 1 cancers among all detected prostate cancers within each study arm."

Plus dla biopsji:

„the ratio of the total number of prostate biopsies performed to the number of men undergoing screening"

Czyli: protokół nazywa dwa mianowniki
mianownik	gdzie w protokole	to jest w analizie
n_screened (przebadani)	§6, §9 + dwa razy w liście	§2.3 tabela 2, §3.1, §3.2, §3.4
pca (wykryte raki)	§6, §9 lista	§2.2, §2.3 tabela 1
n_itt (randomizowani)	nigdzie	cała Part 1 — §2.1, §3.1, §3.2
n_itt nie występuje w protokole ani razu. A to jest mianownik, na którym stoi §2.1 — czyli liczba RR 1.76, którą raport podaje jako wynik główny.

Dwie różne sytuacje, nie jedna
Trzeba to rozdzielić, bo protokół mówi o dwóch rzeczach osobno:

Proporcje (Part 2) — protokół jest jednoznaczny: n_screened i pca. Skrypt robi dokładnie to. ✅ Zgodność.

Risk ratios (Part 1) — protokół §8 mówi tylko „pooled … as risk ratios for comparative outcomes" i nie precyzuje mianownika. Skrypt wybrał n_itt i uzasadnia to w komentarzu:

# RRs always use ITT (randomized / invited) denominators - the randomized
# comparison must not be conditioned on being screened.

To jest metodologicznie słuszna decyzja — i bym jej bronił. Dzieląc przez przebadanych, warunkujesz na zmiennej po randomizacji (bycie przebadanym to decyzja pacjenta: 82% w ERSPC, 35.5% w CAP), więc tracisz ochronę randomizacji — dokładnie ten sam problem, co przy gg1/pca. ITT jest jedynym sposobem, by RR pozostał kauzalny.

Ale to decyzja interpretacyjna, której protokół nie przewidział, a §9 wskazuje wręcz w drugą stronę. Musi trafić do manuskryptu jako zadeklarowane doprecyzowanie. Obecnie REPORT.md §4 deklaruje jedno odchylenie od protokołu (korekta 0.5 przy proporcjach) i o tym nie wspomina.

Napięcie, na które musisz zwrócić uwagę
I tu jest sedno twojego pytania:

mianownik zgodny z protokołem   (n_screened)  →  produkuje artefakt "PSA 0.60%"
mianownik NIEzgodny z protokołem (n_itt)      →  jest metodologicznie poprawny
Protokół każe policzyć „detection rate per screened men". Skrypt to robi (§2.3 tabela 2, §3.1, §3.2) — i wychodzą liczby bez interpretacji, bo mianownik zlicza tylko przebadanych, a licznik zbiera raki z całego ramienia przez 9–10 lat. Raport sam to przyznaje.

Czyli analiza przewidziana protokołem jako pierwotna jest tą zepsutą. To nie wina wykonania — to wada protokołu, napisanego przed zobaczeniem danych.

Który mianownik jest „ten właściwy" — praktycznie
Protokół stawia oba jako współrzędne („among screened men or among detected prostate cancers" — to „or" jest niedoprecyzowaniem po stronie protokołu). Skoro musisz wybrać, dane wybierają za ciebie:

n_screened	pca	n_itt
w protokole	✅ tak	✅ tak	❌ nie
kauzalny	❌ nie (warunkowanie)	❌ nie (warunkowanie)	✅ tak
porównywalny między badaniami	❌ nie (follow-up 9 lat vs 1 runda)	✅ tak (normalizacja wewnątrz badania)	❌ nie
τ² / pooling	🔴 5.90 → artefakt	✅ ≤0.58 → OK	✅ 0.084
Rekomendacja, którą bym wpisał do manuskryptu:

pca jako pierwotny dla tezy pracy (gradient 65.7% → 16.4%). Jest w protokole, τ² małe, odporny na różnice w follow-upie — i raport sam to rozpoznaje w §5 punkcie 4. Opisywać wyłącznie jako skład histologiczny, nie przyczynowo.
n_itt jako pierwotny dla RR (efekt zaproszenia). Zadeklarować jako doprecyzowanie protokołu, z uzasadnieniem ITT.
n_screened jako czysto opisowy, bez pooled estimates tam, gdzie τ² jest duże (grupy PSA i no-screening). Protokół tego wymaga, więc trzeba to policzyć i pokazać — ale jako tabelę arm-level z zakresami, nie jako pooled proportion.
Kluczowe: to jest odejście od protokołu w kwestii hierarchii endpointów, więc wymaga jawnej deklaracji („post-hoc clarification of the primary endpoint"). PRISMA 2020 (pozycja 27) tego wprost wymaga, a protokół §8 sam odwołuje się do PRISMA.

Przy okazji — jeszcze jedna nieścisłość w REPORT.md
Skoro zaglądamy do protokołu, sprawdziłem też to zdanie z §1 raportu:

„two contrast families, which are never pooled into a single interpreted estimate (protocol §8/§9)"

Protokół nie zawiera podziału na rodziny A/B. §8 i §9 mówią wyłącznie:

„Studies will be grouped according to the screening strategy applied, including digital rectal examination (DRE), PSA-based screening, PSA combined with DRE, PSA followed by MRI, MRI-based screening strategies, and other combined or biomarker-based approaches."

To jest grupowanie po strategii (czyli Part 2). Podział „zaproszenie vs brak screeningu" kontra „strategia diagnostyczna przy obu ramionach przebadanych" to dodatek post-hoc.

Dodatek bardzo dobry — bez niego cała analiza byłaby bezsensowna, bo mieszałaby dwa niekompatybilne pytania. Ale przypisanie go protokołowi jest nieprawdziwe. Do przepisania na: „post-hoc stratification by contrast family, introduced because the protocol's strategy grouping does not distinguish invitation trials from diagnostic-strategy trials".

To dopisuję do listy poprawek jako punkt 19.