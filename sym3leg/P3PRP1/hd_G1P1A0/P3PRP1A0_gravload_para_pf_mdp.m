% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRP1G1P1A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [11x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRP1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:42
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRP1G1P1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(11,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1G1P1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1G1P1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRP1G1P1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1G1P1A0_gravload_para_pf_mdp: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1G1P1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1G1P1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [11 1]), ...
  'P3PRP1G1P1A0_gravload_para_pf_mdp: MDP has to be [11x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:42:51
% EndTime: 2019-05-03 14:42:54
% DurationCPUTime: 2.97s
% Computational Cost: add. (1413->248), mult. (2783->399), div. (72->3), fcn. (1282->14), ass. (0->177)
t2741 = -MDP(4) + MDP(6);
t2662 = (pkin(2) ^ 2);
t2733 = t2662 + 1;
t2659 = koppelP(3,1);
t2656 = koppelP(3,2);
t2718 = qJ(3,3) * t2656;
t2610 = pkin(2) * t2659 - t2718;
t2717 = qJ(3,3) * t2659;
t2611 = pkin(2) * t2656 + t2717;
t2642 = legFrame(3,3);
t2628 = sin(t2642);
t2631 = cos(t2642);
t2652 = xP(3);
t2634 = sin(t2652);
t2635 = cos(t2652);
t2737 = (t2610 * t2634 + t2611 * t2635) * t2631 - (t2610 * t2635 - t2611 * t2634) * t2628;
t2660 = koppelP(2,1);
t2657 = koppelP(2,2);
t2723 = qJ(3,2) * t2657;
t2612 = pkin(2) * t2660 - t2723;
t2722 = qJ(3,2) * t2660;
t2613 = pkin(2) * t2657 + t2722;
t2643 = legFrame(2,3);
t2629 = sin(t2643);
t2632 = cos(t2643);
t2736 = (t2612 * t2634 + t2613 * t2635) * t2632 - (t2612 * t2635 - t2613 * t2634) * t2629;
t2661 = koppelP(1,1);
t2658 = koppelP(1,2);
t2728 = qJ(3,1) * t2658;
t2614 = pkin(2) * t2661 - t2728;
t2727 = qJ(3,1) * t2661;
t2615 = pkin(2) * t2658 + t2727;
t2644 = legFrame(1,3);
t2630 = sin(t2644);
t2633 = cos(t2644);
t2735 = (t2614 * t2634 + t2615 * t2635) * t2633 - (t2614 * t2635 - t2615 * t2634) * t2630;
t2734 = pkin(2) * g(2);
t2732 = MDP(3) + MDP(5);
t2731 = qJ(3,1) * t2630;
t2730 = qJ(3,1) * t2633;
t2650 = cos(qJ(2,1));
t2729 = qJ(3,1) * t2650;
t2726 = qJ(3,2) * t2629;
t2725 = qJ(3,2) * t2632;
t2649 = cos(qJ(2,2));
t2724 = qJ(3,2) * t2649;
t2721 = qJ(3,3) * t2628;
t2720 = qJ(3,3) * t2631;
t2648 = cos(qJ(2,3));
t2719 = qJ(3,3) * t2648;
t2645 = sin(qJ(2,3));
t2687 = 0.2e1 * t2719;
t2539 = -t2737 * t2645 + ((t2634 * t2659 + t2635 * t2656) * t2631 - t2628 * (-t2634 * t2656 + t2635 * t2659)) * t2687;
t2653 = qJ(3,3) ^ 2;
t2622 = -t2653 + t2733;
t2639 = t2648 ^ 2;
t2578 = 0.1e1 / (pkin(2) * t2645 * t2687 + t2622 * t2639 - t2653 - t2733);
t2716 = t2539 * t2578;
t2646 = sin(qJ(2,2));
t2689 = 0.2e1 * t2724;
t2540 = -t2736 * t2646 + ((t2634 * t2660 + t2635 * t2657) * t2632 - t2629 * (-t2634 * t2657 + t2635 * t2660)) * t2689;
t2654 = qJ(3,2) ^ 2;
t2623 = -t2654 + t2733;
t2640 = t2649 ^ 2;
t2579 = 0.1e1 / (pkin(2) * t2646 * t2689 + t2623 * t2640 - t2654 - t2733);
t2715 = t2540 * t2579;
t2647 = sin(qJ(2,1));
t2691 = 0.2e1 * t2729;
t2541 = -t2735 * t2647 + ((t2634 * t2661 + t2635 * t2658) * t2633 - t2630 * (-t2634 * t2658 + t2635 * t2661)) * t2691;
t2655 = qJ(3,1) ^ 2;
t2624 = -t2655 + t2733;
t2641 = t2650 ^ 2;
t2580 = 0.1e1 / (pkin(2) * t2647 * t2691 + t2624 * t2641 - t2655 - t2733);
t2714 = t2541 * t2580;
t2688 = -0.2e1 * t2719;
t2569 = t2628 * t2688 + t2645 * (pkin(2) * t2628 - t2720);
t2713 = t2569 * t2578;
t2570 = t2631 * t2688 + t2645 * (pkin(2) * t2631 + t2721);
t2712 = t2570 * t2578;
t2690 = -0.2e1 * t2724;
t2571 = t2629 * t2690 + t2646 * (pkin(2) * t2629 - t2725);
t2711 = t2571 * t2579;
t2572 = t2632 * t2690 + t2646 * (pkin(2) * t2632 + t2726);
t2710 = t2572 * t2579;
t2692 = -0.2e1 * t2729;
t2573 = t2630 * t2692 + t2647 * (pkin(2) * t2630 - t2730);
t2709 = t2573 * t2580;
t2574 = t2633 * t2692 + t2647 * (pkin(2) * t2633 + t2731);
t2708 = t2574 * t2580;
t2590 = g(1) * t2628 - g(2) * t2631;
t2704 = t2578 * t2590;
t2591 = g(1) * t2629 - g(2) * t2632;
t2703 = t2579 * t2591;
t2592 = g(1) * t2630 - g(2) * t2633;
t2702 = t2580 * t2592;
t2701 = t2645 * t2648;
t2700 = t2646 * t2649;
t2699 = t2647 * t2650;
t2698 = t2653 * t2656;
t2697 = t2654 * t2657;
t2696 = t2655 * t2658;
t2695 = t2733 * t2656;
t2694 = t2733 * t2657;
t2693 = t2733 * t2658;
t2686 = pkin(2) * t2727;
t2685 = pkin(2) * t2722;
t2684 = pkin(2) * t2717;
t2683 = pkin(2) * t2731;
t2682 = pkin(2) * t2726;
t2681 = pkin(2) * t2721;
t2680 = pkin(2) * t2720;
t2679 = pkin(2) * t2725;
t2678 = pkin(2) * t2730;
t2625 = pkin(2) * t2718;
t2584 = t2622 * t2659 - 0.2e1 * t2625;
t2585 = 0.2e1 * t2684 + t2695 - t2698;
t2551 = t2584 * t2635 - t2585 * t2634;
t2552 = t2584 * t2634 + t2585 * t2635;
t2626 = pkin(2) * t2723;
t2586 = t2623 * t2660 - 0.2e1 * t2626;
t2587 = 0.2e1 * t2685 + t2694 - t2697;
t2553 = t2586 * t2635 - t2587 * t2634;
t2554 = t2586 * t2634 + t2587 * t2635;
t2627 = pkin(2) * t2728;
t2588 = t2624 * t2661 - 0.2e1 * t2627;
t2589 = 0.2e1 * t2686 + t2693 - t2696;
t2555 = t2588 * t2635 - t2589 * t2634;
t2556 = t2588 * t2634 + t2589 * t2635;
t2604 = t2659 * t2662 - t2625 + t2659;
t2605 = t2684 + t2695;
t2606 = t2660 * t2662 - t2626 + t2660;
t2607 = t2685 + t2694;
t2608 = t2661 * t2662 - t2627 + t2661;
t2609 = t2686 + t2693;
t2677 = ((t2555 * t2633 + t2556 * t2630) * t2641 + (-t2555 * t2630 + t2556 * t2633) * t2699 + (-t2608 * t2635 + t2609 * t2634) * t2633 - (t2608 * t2634 + t2609 * t2635) * t2630) * t2702 + ((t2553 * t2632 + t2554 * t2629) * t2640 + (-t2553 * t2629 + t2554 * t2632) * t2700 + (-t2606 * t2635 + t2607 * t2634) * t2632 - (t2606 * t2634 + t2607 * t2635) * t2629) * t2703 + ((t2551 * t2631 + t2552 * t2628) * t2639 + (-t2551 * t2628 + t2552 * t2631) * t2701 + (-t2604 * t2635 + t2605 * t2634) * t2631 - (t2604 * t2634 + t2605 * t2635) * t2628) * t2704;
t2581 = t2622 * t2631 + 0.2e1 * t2681;
t2582 = t2623 * t2632 + 0.2e1 * t2682;
t2583 = t2624 * t2633 + 0.2e1 * t2683;
t2665 = -t2622 * t2628 + 0.2e1 * t2680;
t2668 = -t2623 * t2629 + 0.2e1 * t2679;
t2671 = -t2624 * t2630 + 0.2e1 * t2678;
t2676 = (t2583 * t2641 - t2633 * t2662 + t2671 * t2699 - t2633 - t2683) * t2702 + (t2582 * t2640 - t2632 * t2662 + t2668 * t2700 - t2632 - t2682) * t2703 + (t2581 * t2639 - t2631 * t2662 + t2665 * t2701 - t2631 - t2681) * t2704;
t2675 = (-t2583 * t2699 + t2662 * t2630 + t2671 * t2641 + t2630 - t2678) * t2702 + (-t2582 * t2700 + t2662 * t2629 + t2668 * t2640 + t2629 - t2679) * t2703 + (-t2581 * t2701 + t2662 * t2628 + t2665 * t2639 + t2628 - t2680) * t2704;
t2593 = g(1) * t2631 + g(2) * t2628;
t2560 = t2590 * t2648 + t2593 * t2645;
t2557 = t2590 * t2645 - t2593 * t2648;
t2594 = g(1) * t2632 + g(2) * t2629;
t2561 = t2591 * t2649 + t2594 * t2646;
t2558 = t2591 * t2646 - t2594 * t2649;
t2595 = g(1) * t2633 + g(2) * t2630;
t2562 = t2592 * t2650 + t2595 * t2647;
t2559 = t2592 * t2647 - t2595 * t2650;
t2670 = -t2630 * t2655 - t2678;
t2669 = t2633 * t2655 - t2683;
t2667 = -t2629 * t2654 - t2679;
t2666 = t2632 * t2654 - t2682;
t2664 = -t2628 * t2653 - t2680;
t2663 = t2631 * t2653 - t2681;
t2651 = pkin(2) * g(1);
t2621 = g(1) * qJ(3,1) + t2734;
t2620 = -g(2) * qJ(3,1) + t2651;
t2619 = g(1) * qJ(3,2) + t2734;
t2618 = -g(2) * qJ(3,2) + t2651;
t2617 = g(1) * qJ(3,3) + t2734;
t2616 = -g(2) * qJ(3,3) + t2651;
t2603 = t2655 * t2661 + t2627 + t2661;
t2602 = -t2658 + t2686 - t2696;
t2601 = t2654 * t2660 + t2626 + t2660;
t2600 = -t2657 + t2685 - t2697;
t2599 = t2653 * t2659 + t2625 + t2659;
t2598 = -t2656 + t2684 - t2698;
t2597 = g(1) * t2635 + g(2) * t2634;
t2596 = g(1) * t2634 - g(2) * t2635;
t2550 = (t2620 * t2647 - t2621 * t2650) * t2633 + (t2620 * t2650 + t2621 * t2647) * t2630;
t2549 = (t2618 * t2646 - t2619 * t2649) * t2632 + (t2618 * t2649 + t2619 * t2646) * t2629;
t2548 = (t2616 * t2645 - t2617 * t2648) * t2631 + (t2616 * t2648 + t2617 * t2645) * t2628;
t1 = [t2675 * MDP(1) + ((t2574 * t2550 - (t2670 * t2650 + t2647 * (-t2633 - t2669)) * t2562) * t2580 + (t2572 * t2549 - (t2667 * t2649 + t2646 * (-t2632 - t2666)) * t2561) * t2579 + (t2570 * t2548 - (t2664 * t2648 + t2645 * (-t2631 - t2663)) * t2560) * t2578 + t2675) * MDP(7) + (-t2596 * t2634 - t2597 * t2635) * MDP(11) + t2732 * (t2560 * t2712 + t2561 * t2710 + t2562 * t2708) + t2741 * (t2557 * t2712 + t2558 * t2710 + t2559 * t2708); t2676 * MDP(1) + ((t2573 * t2550 - (t2669 * t2650 - t2647 * (t2630 - t2670)) * t2562) * t2580 + (t2571 * t2549 - (t2666 * t2649 - t2646 * (t2629 - t2667)) * t2561) * t2579 + (t2569 * t2548 - (t2663 * t2648 - t2645 * (t2628 - t2664)) * t2560) * t2578 + t2676) * MDP(7) + (t2596 * t2635 - t2597 * t2634) * MDP(11) + t2732 * (t2560 * t2713 + t2561 * t2711 + t2562 * t2709) + t2741 * (t2557 * t2713 + t2558 * t2711 + t2559 * t2709); t2677 * MDP(1) + ((t2541 * t2550 - (((-t2602 * t2635 + t2603 * t2634) * t2633 - (t2602 * t2634 + t2603 * t2635) * t2630) * t2647 + t2735 * t2729) * t2562) * t2580 + (t2540 * t2549 - (((-t2600 * t2635 + t2601 * t2634) * t2632 - (t2600 * t2634 + t2601 * t2635) * t2629) * t2646 + t2736 * t2724) * t2561) * t2579 + (t2539 * t2548 - (((-t2598 * t2635 + t2599 * t2634) * t2631 - (t2598 * t2634 + t2599 * t2635) * t2628) * t2645 + t2737 * t2719) * t2560) * t2578 + t2677) * MDP(7) + t2596 * MDP(9) + t2597 * MDP(10) + t2732 * (t2560 * t2716 + t2561 * t2715 + t2562 * t2714) + t2741 * (t2557 * t2716 + t2558 * t2715 + t2559 * t2714);];
taugX  = t1;
