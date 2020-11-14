% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR8V1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR8V1G1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 16:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V1G1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 16:50:51
% EndTime: 2020-08-06 16:50:52
% DurationCPUTime: 1.31s
% Computational Cost: add. (661->183), mult. (1635->366), div. (84->7), fcn. (1908->22), ass. (0->139)
t2618 = cos(qJ(3,1));
t2619 = cos(qJ(2,1));
t2613 = sin(qJ(2,1));
t2645 = t2613 * t2618;
t2582 = pkin(2) * t2645 - pkin(5) * t2619;
t2602 = sin(pkin(3));
t2604 = cos(pkin(3));
t2612 = sin(qJ(3,1));
t2658 = t2604 * t2612;
t2627 = 0.1e1 / (pkin(2) * t2658 + t2582 * t2602);
t2685 = t2627 / t2618;
t2616 = cos(qJ(3,2));
t2617 = cos(qJ(2,2));
t2611 = sin(qJ(2,2));
t2648 = t2611 * t2616;
t2581 = pkin(2) * t2648 - pkin(5) * t2617;
t2610 = sin(qJ(3,2));
t2660 = t2604 * t2610;
t2628 = 0.1e1 / (pkin(2) * t2660 + t2581 * t2602);
t2684 = t2628 / t2616;
t2614 = cos(qJ(3,3));
t2615 = cos(qJ(2,3));
t2609 = sin(qJ(2,3));
t2651 = t2609 * t2614;
t2580 = pkin(2) * t2651 - pkin(5) * t2615;
t2608 = sin(qJ(3,3));
t2662 = t2604 * t2608;
t2629 = 0.1e1 / (pkin(2) * t2662 + t2580 * t2602);
t2683 = t2629 / t2614;
t2682 = MDP(1) * g(3);
t2681 = pkin(2) * t2614;
t2680 = pkin(2) * t2616;
t2679 = pkin(2) * t2618;
t2678 = g(3) * t2602;
t2601 = sin(pkin(6));
t2603 = cos(pkin(6));
t2586 = g(1) * t2601 - g(2) * t2603;
t2587 = g(1) * t2603 + g(2) * t2601;
t2605 = legFrame(3,3);
t2592 = sin(t2605);
t2595 = cos(t2605);
t2529 = (-t2678 + (t2586 * t2595 + t2587 * t2592) * t2604) * t2615 + (-t2586 * t2592 + t2587 * t2595) * t2609;
t2677 = t2529 * t2629;
t2606 = legFrame(2,3);
t2593 = sin(t2606);
t2596 = cos(t2606);
t2530 = (-t2678 + (t2586 * t2596 + t2587 * t2593) * t2604) * t2617 + (-t2586 * t2593 + t2587 * t2596) * t2611;
t2676 = t2530 * t2628;
t2607 = legFrame(1,3);
t2594 = sin(t2607);
t2597 = cos(t2607);
t2531 = (-t2678 + (t2586 * t2597 + t2587 * t2594) * t2604) * t2619 + (-t2586 * t2594 + t2587 * t2597) * t2613;
t2675 = t2531 * t2627;
t2668 = t2602 * t2608;
t2667 = t2602 * t2610;
t2666 = t2602 * t2612;
t2665 = t2602 * t2614;
t2664 = t2602 * t2616;
t2663 = t2602 * t2618;
t2661 = t2604 * t2609;
t2659 = t2604 * t2611;
t2657 = t2604 * t2613;
t2656 = t2604 * t2615;
t2655 = t2604 * t2617;
t2654 = t2604 * t2619;
t2653 = t2608 * t2609;
t2652 = t2608 * t2615;
t2650 = t2610 * t2611;
t2649 = t2610 * t2617;
t2647 = t2612 * t2613;
t2646 = t2612 * t2619;
t2644 = t2614 * t2615;
t2643 = t2616 * t2617;
t2642 = t2618 * t2619;
t2556 = -t2592 * t2601 + t2595 * t2603;
t2559 = t2592 * t2603 + t2595 * t2601;
t2641 = (-(t2556 * t2656 - t2559 * t2609) * t2681 - pkin(5) * (t2556 * t2661 + t2559 * t2615)) * t2683;
t2640 = (-(t2556 * t2609 + t2559 * t2656) * t2681 - (-t2556 * t2615 + t2559 * t2661) * pkin(5)) * t2683;
t2557 = -t2593 * t2601 + t2596 * t2603;
t2560 = t2593 * t2603 + t2596 * t2601;
t2639 = (-(t2557 * t2655 - t2560 * t2611) * t2680 - pkin(5) * (t2557 * t2659 + t2560 * t2617)) * t2684;
t2638 = (-(t2557 * t2611 + t2560 * t2655) * t2680 - (-t2557 * t2617 + t2560 * t2659) * pkin(5)) * t2684;
t2558 = -t2594 * t2601 + t2597 * t2603;
t2561 = t2594 * t2603 + t2597 * t2601;
t2637 = (-(t2558 * t2654 - t2561 * t2613) * t2679 - pkin(5) * (t2558 * t2657 + t2561 * t2619)) * t2685;
t2636 = (-(t2558 * t2613 + t2561 * t2654) * t2679 - (-t2558 * t2619 + t2561 * t2657) * pkin(5)) * t2685;
t2574 = g(1) * t2592 - g(2) * t2595;
t2577 = g(1) * t2595 + g(2) * t2592;
t2635 = (t2577 * (t2601 * t2656 + t2603 * t2609) + t2574 * (-t2601 * t2609 + t2603 * t2656) - t2615 * t2678) * t2683;
t2563 = t2601 * t2615 + t2603 * t2661;
t2566 = -t2601 * t2661 + t2603 * t2615;
t2634 = (-t2563 * t2574 + t2566 * t2577 + t2609 * t2678) * t2683;
t2575 = g(1) * t2593 - g(2) * t2596;
t2578 = g(1) * t2596 + g(2) * t2593;
t2633 = (t2578 * (t2601 * t2655 + t2603 * t2611) + t2575 * (-t2601 * t2611 + t2603 * t2655) - t2617 * t2678) * t2684;
t2564 = t2601 * t2617 + t2603 * t2659;
t2567 = -t2601 * t2659 + t2603 * t2617;
t2632 = (-t2564 * t2575 + t2567 * t2578 + t2611 * t2678) * t2684;
t2576 = g(1) * t2594 - g(2) * t2597;
t2579 = g(1) * t2597 + g(2) * t2594;
t2631 = (t2579 * (t2601 * t2654 + t2603 * t2613) + t2576 * (-t2601 * t2613 + t2603 * t2654) - t2619 * t2678) * t2685;
t2562 = t2601 * t2657 - t2603 * t2619;
t2565 = t2601 * t2619 + t2603 * t2657;
t2630 = (-t2562 * t2579 - t2565 * t2576 + t2613 * t2678) * t2685;
t2626 = t2529 * t2608 * t2683;
t2625 = t2530 * t2610 * t2684;
t2624 = t2531 * t2612 * t2685;
t2623 = pkin(2) * t2668 - t2580 * t2604;
t2622 = pkin(2) * t2667 - t2581 * t2604;
t2621 = pkin(2) * t2666 - t2582 * t2604;
t2620 = 0.1e1 / pkin(2);
t2585 = pkin(2) * t2642 + pkin(5) * t2613;
t2584 = pkin(2) * t2643 + pkin(5) * t2611;
t2583 = pkin(2) * t2644 + pkin(5) * t2609;
t2573 = t2604 * t2645 - t2666;
t2572 = t2604 * t2647 + t2663;
t2571 = t2604 * t2648 - t2667;
t2570 = t2604 * t2651 - t2668;
t2569 = t2604 * t2650 + t2664;
t2568 = t2604 * t2653 + t2665;
t2549 = t2585 * t2601 - t2621 * t2603;
t2548 = t2584 * t2601 - t2622 * t2603;
t2547 = t2583 * t2601 - t2623 * t2603;
t2546 = t2585 * t2603 + t2621 * t2601;
t2545 = t2584 * t2603 + t2622 * t2601;
t2544 = t2583 * t2603 + t2623 * t2601;
t2537 = (-t2562 * t2597 - t2565 * t2594) * t2612 - t2561 * t2663;
t2536 = (-t2564 * t2593 + t2567 * t2596) * t2610 - t2560 * t2664;
t2535 = (-t2563 * t2592 + t2566 * t2595) * t2608 - t2559 * t2665;
t2534 = (t2562 * t2594 - t2565 * t2597) * t2612 - t2558 * t2663;
t2533 = (-t2564 * t2596 - t2567 * t2593) * t2610 - t2557 * t2664;
t2532 = (-t2563 * t2595 - t2566 * t2592) * t2608 - t2556 * t2665;
t2522 = (-t2573 * t2601 + t2603 * t2642) * t2579 - t2576 * (t2573 * t2603 + t2601 * t2642) + g(3) * (t2602 * t2645 + t2658);
t2521 = (-t2571 * t2601 + t2603 * t2643) * t2578 - t2575 * (t2571 * t2603 + t2601 * t2643) + g(3) * (t2602 * t2648 + t2660);
t2520 = (-t2570 * t2601 + t2603 * t2644) * t2577 - t2574 * (t2570 * t2603 + t2601 * t2644) + g(3) * (t2602 * t2651 + t2662);
t2519 = t2579 * (-t2572 * t2601 + t2603 * t2646) - (t2572 * t2603 + t2601 * t2646) * t2576 - g(3) * (-t2602 * t2647 + t2604 * t2618);
t2518 = t2578 * (-t2569 * t2601 + t2603 * t2649) - (t2569 * t2603 + t2601 * t2649) * t2575 - g(3) * (-t2602 * t2650 + t2604 * t2616);
t2517 = t2577 * (-t2568 * t2601 + t2603 * t2652) - (t2568 * t2603 + t2601 * t2652) * t2574 - g(3) * (-t2602 * t2653 + t2604 * t2614);
t1 = [(-(t2546 * t2597 - t2549 * t2594) * t2627 - (t2545 * t2596 - t2548 * t2593) * t2628 - (t2544 * t2595 - t2547 * t2592) * t2629) * t2682 + (t2532 * t2635 + t2533 * t2633 + t2534 * t2631) * MDP(3) + (t2532 * t2634 + t2533 * t2632 + t2534 * t2630) * MDP(4) + (t2532 * t2677 + t2533 * t2676 + t2534 * t2675 + (t2517 * t2641 + t2518 * t2639 + t2519 * t2637) * t2620) * MDP(10) + (-t2532 * t2626 - t2533 * t2625 - t2534 * t2624 + (t2520 * t2641 + t2521 * t2639 + t2522 * t2637) * t2620) * MDP(11) - g(1) * MDP(12); (-(t2546 * t2594 + t2549 * t2597) * t2627 - (t2545 * t2593 + t2548 * t2596) * t2628 - (t2544 * t2592 + t2547 * t2595) * t2629) * t2682 + (t2535 * t2635 + t2536 * t2633 + t2537 * t2631) * MDP(3) + (t2535 * t2634 + t2536 * t2632 + t2537 * t2630) * MDP(4) + (t2535 * t2677 + t2536 * t2676 + t2537 * t2675 + (t2517 * t2640 + t2518 * t2638 + t2519 * t2636) * t2620) * MDP(10) + (-t2535 * t2626 - t2536 * t2625 - t2537 * t2624 + (t2520 * t2640 + t2521 * t2638 + t2522 * t2636) * t2620) * MDP(11) - g(2) * MDP(12); (-(3 * MDP(1)) - MDP(12)) * g(3);];
taugX  = t1;
