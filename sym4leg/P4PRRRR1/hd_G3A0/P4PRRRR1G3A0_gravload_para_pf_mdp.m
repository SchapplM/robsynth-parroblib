% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P4PRRRR1G3P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P4PRRRR1G3P1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [4x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 19:06
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR1G3P1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 19:06:23
% EndTime: 2020-03-02 19:06:25
% DurationCPUTime: 2.58s
% Computational Cost: add. (413->153), mult. (829->297), div. (220->13), fcn. (932->26), ass. (0->131)
t2564 = legFrame(4,2);
t2538 = sin(t2564);
t2542 = cos(t2564);
t2528 = g(1) * t2538 + g(2) * t2542;
t2532 = g(1) * t2542 - g(2) * t2538;
t2561 = sin(qJ(2,4));
t2563 = cos(qJ(2,4));
t2506 = t2528 * t2563 - t2532 * t2561;
t2548 = 0.1e1 / t2561;
t2637 = t2506 * t2548;
t2565 = legFrame(3,2);
t2539 = sin(t2565);
t2543 = cos(t2565);
t2529 = g(1) * t2539 + g(2) * t2543;
t2533 = g(1) * t2543 - g(2) * t2539;
t2569 = sin(qJ(2,3));
t2575 = cos(qJ(2,3));
t2509 = t2529 * t2575 - t2533 * t2569;
t2551 = 0.1e1 / t2569;
t2636 = t2509 * t2551;
t2566 = legFrame(2,2);
t2540 = sin(t2566);
t2544 = cos(t2566);
t2530 = g(1) * t2540 + g(2) * t2544;
t2534 = g(1) * t2544 - g(2) * t2540;
t2571 = sin(qJ(2,2));
t2577 = cos(qJ(2,2));
t2512 = t2530 * t2577 - t2534 * t2571;
t2552 = 0.1e1 / t2571;
t2635 = t2512 * t2552;
t2567 = legFrame(1,2);
t2541 = sin(t2567);
t2545 = cos(t2567);
t2531 = g(1) * t2541 + g(2) * t2545;
t2535 = g(1) * t2545 - g(2) * t2541;
t2573 = sin(qJ(2,1));
t2579 = cos(qJ(2,1));
t2515 = t2531 * t2579 - t2535 * t2573;
t2553 = 0.1e1 / t2573;
t2634 = t2515 * t2553;
t2633 = t2528 * t2548;
t2632 = t2529 * t2551;
t2631 = t2530 * t2552;
t2630 = t2531 * t2553;
t2562 = cos(qJ(3,4));
t2549 = 0.1e1 / t2562;
t2629 = t2548 * t2549;
t2574 = cos(qJ(3,3));
t2554 = 0.1e1 / t2574;
t2628 = t2551 * t2554;
t2576 = cos(qJ(3,2));
t2556 = 0.1e1 / t2576;
t2627 = t2552 * t2556;
t2578 = cos(qJ(3,1));
t2558 = 0.1e1 / t2578;
t2626 = t2553 * t2558;
t2580 = xP(4);
t2546 = sin(t2580);
t2547 = cos(t2580);
t2581 = koppelP(4,2);
t2585 = koppelP(4,1);
t2596 = -t2546 * t2581 + t2547 * t2585;
t2597 = t2546 * t2585 + t2547 * t2581;
t2500 = t2596 * t2538 + t2597 * t2542;
t2625 = t2500 * t2629;
t2582 = koppelP(3,2);
t2586 = koppelP(3,1);
t2594 = -t2546 * t2582 + t2547 * t2586;
t2595 = t2546 * t2586 + t2547 * t2582;
t2501 = t2594 * t2539 + t2595 * t2543;
t2624 = t2501 * t2628;
t2583 = koppelP(2,2);
t2587 = koppelP(2,1);
t2592 = -t2546 * t2583 + t2547 * t2587;
t2593 = t2546 * t2587 + t2547 * t2583;
t2502 = t2592 * t2540 + t2593 * t2544;
t2623 = t2502 * t2627;
t2584 = koppelP(1,2);
t2588 = koppelP(1,1);
t2590 = -t2546 * t2584 + t2547 * t2588;
t2591 = t2546 * t2588 + t2547 * t2584;
t2503 = t2590 * t2541 + t2591 * t2545;
t2622 = t2503 * t2626;
t2621 = t2538 * t2629;
t2620 = t2539 * t2628;
t2619 = t2540 * t2627;
t2618 = t2541 * t2626;
t2617 = t2542 * t2629;
t2616 = t2543 * t2628;
t2615 = t2544 * t2627;
t2614 = t2545 * t2626;
t2560 = sin(qJ(3,4));
t2613 = t2560 * t2629;
t2612 = t2548 / t2562 ^ 2 * t2563;
t2568 = sin(qJ(3,3));
t2611 = t2568 * t2628;
t2610 = t2551 / t2574 ^ 2 * t2575;
t2570 = sin(qJ(3,2));
t2609 = t2570 * t2627;
t2608 = t2552 / t2576 ^ 2 * t2577;
t2572 = sin(qJ(3,1));
t2607 = t2572 * t2626;
t2606 = t2553 / t2578 ^ 2 * t2579;
t2605 = t2506 * t2613;
t2604 = t2509 * t2611;
t2603 = t2512 * t2609;
t2602 = t2515 * t2607;
t2601 = t2560 * t2612;
t2600 = t2568 * t2610;
t2599 = t2570 * t2608;
t2598 = t2572 * t2606;
t2589 = 0.1e1 / pkin(2);
t2537 = g(1) * t2547 + g(2) * t2546;
t2536 = g(1) * t2546 - g(2) * t2547;
t2527 = t2541 * t2573 + t2545 * t2579;
t2526 = -t2541 * t2579 + t2545 * t2573;
t2525 = t2540 * t2571 + t2544 * t2577;
t2524 = -t2540 * t2577 + t2544 * t2571;
t2523 = t2539 * t2569 + t2543 * t2575;
t2522 = -t2539 * t2575 + t2543 * t2569;
t2521 = t2538 * t2561 + t2542 * t2563;
t2520 = -t2538 * t2563 + t2542 * t2561;
t2519 = (g(1) * t2579 + g(2) * t2573) * t2545 + t2541 * (g(1) * t2573 - g(2) * t2579);
t2518 = (g(1) * t2577 + g(2) * t2571) * t2544 + t2540 * (g(1) * t2571 - g(2) * t2577);
t2517 = (g(1) * t2575 + g(2) * t2569) * t2543 + t2539 * (g(1) * t2569 - g(2) * t2575);
t2516 = (g(1) * t2563 + g(2) * t2561) * t2542 + t2538 * (g(1) * t2561 - g(2) * t2563);
t2514 = t2531 * t2573 + t2535 * t2579;
t2511 = t2530 * t2571 + t2534 * t2577;
t2508 = t2529 * t2569 + t2533 * t2575;
t2505 = t2528 * t2561 + t2532 * t2563;
t1 = [(-t2521 * t2633 - t2523 * t2632 - t2525 * t2631 - t2527 * t2630) * MDP(1) + (-t2536 * t2546 - t2537 * t2547) * MDP(15) + ((t2506 * t2617 + t2509 * t2616 + t2512 * t2615 + t2515 * t2614) * MDP(3) + (-t2505 * t2617 - t2508 * t2616 - t2511 * t2615 - t2514 * t2614) * MDP(4) + (t2542 * t2637 + t2543 * t2636 + t2544 * t2635 + t2545 * t2634) * MDP(10) + (-t2542 * t2605 - t2543 * t2604 - t2544 * t2603 - t2545 * t2602) * MDP(11)) * t2589; (-t2520 * t2633 - t2522 * t2632 - t2524 * t2631 - t2526 * t2630) * MDP(1) + (t2536 * t2547 - t2537 * t2546) * MDP(15) + ((-t2506 * t2621 - t2509 * t2620 - t2512 * t2619 - t2515 * t2618) * MDP(3) + (t2505 * t2621 + t2508 * t2620 + t2511 * t2619 + t2514 * t2618) * MDP(4) + (-t2538 * t2637 - t2539 * t2636 - t2540 * t2635 - t2541 * t2634) * MDP(10) + (t2538 * t2605 + t2539 * t2604 + t2540 * t2603 + t2541 * t2602) * MDP(11)) * t2589; (-t2528 * t2613 - t2529 * t2611 - t2530 * t2609 - t2531 * t2607) * MDP(1) - g(3) * MDP(15) + ((t2506 * t2601 + t2509 * t2600 + t2512 * t2599 + t2515 * t2598) * MDP(3) + (-t2505 * t2601 - t2508 * t2600 - t2511 * t2599 - t2514 * t2598) * MDP(4) + (t2579 * t2602 + t2558 * (-g(3) * t2578 + t2519 * t2572) + t2577 * t2603 + t2556 * (-g(3) * t2576 + t2518 * t2570) + t2575 * t2604 + t2554 * (-g(3) * t2574 + t2517 * t2568) + t2563 * t2605 + t2549 * (-g(3) * t2562 + t2516 * t2560)) * MDP(10) + (-t2572 ^ 2 * t2515 * t2606 + t2558 * (g(3) * t2572 + t2519 * t2578) - t2570 ^ 2 * t2512 * t2608 + t2556 * (g(3) * t2570 + t2518 * t2576) - t2568 ^ 2 * t2509 * t2610 + t2554 * (g(3) * t2568 + t2517 * t2574) - t2560 ^ 2 * t2506 * t2612 + t2549 * (g(3) * t2560 + t2516 * t2562)) * MDP(11)) * t2589; (-(t2526 * t2590 - t2527 * t2591) * t2630 - (t2524 * t2592 - t2525 * t2593) * t2631 - (t2522 * t2594 - t2523 * t2595) * t2632 - (t2520 * t2596 - t2521 * t2597) * t2633) * MDP(1) + t2536 * MDP(13) + t2537 * MDP(14) + ((-t2506 * t2625 - t2509 * t2624 - t2512 * t2623 - t2515 * t2622) * MDP(3) + (t2505 * t2625 + t2508 * t2624 + t2511 * t2623 + t2514 * t2622) * MDP(4) + (-t2500 * t2637 - t2501 * t2636 - t2502 * t2635 - t2503 * t2634) * MDP(10) + (t2500 * t2605 + t2501 * t2604 + t2502 * t2603 + t2503 * t2602) * MDP(11)) * t2589;];
taugX  = t1;
