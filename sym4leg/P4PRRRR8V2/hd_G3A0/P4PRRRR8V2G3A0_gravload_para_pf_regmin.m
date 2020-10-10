% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P4PRRRR8V2G3A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [4x15]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P4PRRRR8V2G3A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_regmin: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:28:28
% EndTime: 2020-08-07 11:28:32
% DurationCPUTime: 4.97s
% Computational Cost: add. (2113->374), mult. (4920->737), div. (148->9), fcn. (5124->30), ass. (0->293)
t2649 = sin(pkin(4));
t2670 = cos(qJ(3,1));
t2665 = sin(qJ(2,1));
t2831 = pkin(3) * t2670 ^ 2;
t2747 = t2665 * t2831;
t2671 = cos(qJ(2,1));
t2672 = pkin(7) + pkin(6);
t2843 = pkin(2) * t2665 - t2671 * t2672;
t2850 = (t2670 * t2843 + t2747) * t2649;
t2668 = cos(qJ(3,2));
t2663 = sin(qJ(2,2));
t2832 = pkin(3) * t2668 ^ 2;
t2748 = t2663 * t2832;
t2669 = cos(qJ(2,2));
t2844 = pkin(2) * t2663 - t2669 * t2672;
t2849 = (t2668 * t2844 + t2748) * t2649;
t2666 = cos(qJ(3,3));
t2661 = sin(qJ(2,3));
t2833 = pkin(3) * t2666 ^ 2;
t2749 = t2661 * t2833;
t2667 = cos(qJ(2,3));
t2845 = pkin(2) * t2661 - t2667 * t2672;
t2848 = (t2666 * t2845 + t2749) * t2649;
t2654 = cos(qJ(3,4));
t2653 = sin(qJ(2,4));
t2834 = pkin(3) * t2654 ^ 2;
t2750 = t2653 * t2834;
t2655 = cos(qJ(2,4));
t2846 = pkin(2) * t2653 - t2655 * t2672;
t2847 = (t2654 * t2846 + t2750) * t2649;
t2651 = cos(pkin(4));
t2826 = t2670 * pkin(3);
t2743 = pkin(2) + t2826;
t2700 = t2651 * t2743;
t2828 = t2666 * pkin(3);
t2745 = pkin(2) + t2828;
t2701 = t2651 * t2745;
t2829 = t2654 * pkin(3);
t2746 = pkin(2) + t2829;
t2702 = t2651 * t2746;
t2827 = t2668 * pkin(3);
t2744 = pkin(2) + t2827;
t2699 = t2744 * t2651;
t2652 = sin(qJ(3,4));
t2838 = pkin(2) * t2652;
t2660 = sin(qJ(3,3));
t2837 = pkin(2) * t2660;
t2662 = sin(qJ(3,2));
t2836 = pkin(2) * t2662;
t2664 = sin(qJ(3,1));
t2835 = pkin(2) * t2664;
t2648 = sin(pkin(8));
t2830 = t2648 * g(3);
t2785 = t2651 * t2652;
t2568 = pkin(3) * t2785 + t2649 * t2846;
t2544 = 0.1e1 / (pkin(2) * t2785 + t2568 * t2654 + t2649 * t2750);
t2796 = t2648 * t2651;
t2602 = -t2649 * g(1) - g(2) * t2796;
t2603 = g(1) * t2796 - t2649 * g(2);
t2656 = legFrame(4,2);
t2631 = sin(t2656);
t2635 = cos(t2656);
t2608 = t2635 * g(1) - t2631 * g(2);
t2650 = cos(pkin(8));
t2790 = t2650 * t2651;
t2625 = g(3) * t2790;
t2825 = ((t2602 * t2631 + t2603 * t2635 + t2625) * t2655 + t2653 * (t2608 * t2650 - t2830)) * t2544;
t2782 = t2651 * t2660;
t2569 = pkin(3) * t2782 + t2649 * t2845;
t2545 = 0.1e1 / (pkin(2) * t2782 + t2569 * t2666 + t2649 * t2749);
t2657 = legFrame(3,2);
t2632 = sin(t2657);
t2636 = cos(t2657);
t2609 = t2636 * g(1) - t2632 * g(2);
t2824 = ((t2602 * t2632 + t2603 * t2636 + t2625) * t2667 + t2661 * (t2609 * t2650 - t2830)) * t2545;
t2780 = t2651 * t2662;
t2570 = pkin(3) * t2780 + t2649 * t2844;
t2546 = 0.1e1 / (pkin(2) * t2780 + t2570 * t2668 + t2649 * t2748);
t2658 = legFrame(2,2);
t2633 = sin(t2658);
t2637 = cos(t2658);
t2610 = t2637 * g(1) - t2633 * g(2);
t2823 = ((t2602 * t2633 + t2603 * t2637 + t2625) * t2669 + t2663 * (t2610 * t2650 - t2830)) * t2546;
t2778 = t2651 * t2664;
t2571 = pkin(3) * t2778 + t2649 * t2843;
t2547 = 0.1e1 / (pkin(2) * t2778 + t2571 * t2670 + t2649 * t2747);
t2659 = legFrame(1,2);
t2634 = sin(t2659);
t2638 = cos(t2659);
t2611 = t2638 * g(1) - t2634 * g(2);
t2822 = ((t2602 * t2634 + t2603 * t2638 + t2625) * t2671 + t2665 * (t2611 * t2650 - t2830)) * t2547;
t2783 = t2651 * t2655;
t2575 = t2648 * t2783 + t2650 * t2653;
t2615 = pkin(2) * t2655 + t2653 * t2672;
t2821 = (t2575 * t2829 + t2615 * t2796 + t2650 * t2846) * t2544;
t2776 = t2651 * t2667;
t2586 = t2648 * t2776 + t2650 * t2661;
t2619 = pkin(2) * t2667 + t2661 * t2672;
t2820 = (t2586 * t2828 + t2619 * t2796 + t2650 * t2845) * t2545;
t2775 = t2651 * t2669;
t2587 = t2648 * t2775 + t2650 * t2663;
t2620 = pkin(2) * t2669 + t2663 * t2672;
t2819 = (t2587 * t2827 + t2620 * t2796 + t2650 * t2844) * t2546;
t2774 = t2651 * t2671;
t2588 = t2648 * t2774 + t2650 * t2665;
t2621 = pkin(2) * t2671 + t2665 * t2672;
t2818 = (t2588 * t2826 + t2621 * t2796 + t2650 * t2843) * t2547;
t2673 = xP(4);
t2642 = sin(t2673);
t2643 = cos(t2673);
t2674 = koppelP(4,2);
t2678 = koppelP(4,1);
t2594 = t2642 * t2678 + t2643 * t2674;
t2598 = -t2642 * t2674 + t2643 * t2678;
t2553 = t2594 * t2635 + t2631 * t2598;
t2817 = t2544 * t2553;
t2784 = t2651 * t2653;
t2789 = t2650 * t2655;
t2572 = t2648 * t2784 - t2789;
t2769 = t2654 * t2649;
t2560 = t2652 * t2572 + t2648 * t2769;
t2816 = t2544 * t2560;
t2604 = t2631 * g(1) + t2635 * g(2);
t2815 = t2544 * t2604;
t2814 = t2544 * t2631;
t2813 = t2544 * t2635;
t2675 = koppelP(3,2);
t2679 = koppelP(3,1);
t2595 = t2642 * t2679 + t2643 * t2675;
t2599 = -t2642 * t2675 + t2643 * t2679;
t2554 = t2595 * t2636 + t2632 * t2599;
t2812 = t2545 * t2554;
t2781 = t2651 * t2661;
t2788 = t2650 * t2667;
t2577 = t2648 * t2781 - t2788;
t2756 = t2666 * t2649;
t2563 = t2660 * t2577 + t2648 * t2756;
t2811 = t2545 * t2563;
t2605 = t2632 * g(1) + t2636 * g(2);
t2810 = t2545 * t2605;
t2809 = t2545 * t2632;
t2808 = t2545 * t2636;
t2676 = koppelP(2,2);
t2680 = koppelP(2,1);
t2596 = t2642 * t2680 + t2643 * t2676;
t2600 = -t2642 * t2676 + t2643 * t2680;
t2555 = t2596 * t2637 + t2633 * t2600;
t2807 = t2546 * t2555;
t2779 = t2651 * t2663;
t2787 = t2650 * t2669;
t2578 = t2648 * t2779 - t2787;
t2754 = t2668 * t2649;
t2564 = t2662 * t2578 + t2648 * t2754;
t2806 = t2546 * t2564;
t2606 = t2633 * g(1) + t2637 * g(2);
t2805 = t2546 * t2606;
t2804 = t2546 * t2633;
t2803 = t2546 * t2637;
t2677 = koppelP(1,2);
t2681 = koppelP(1,1);
t2597 = t2642 * t2681 + t2643 * t2677;
t2601 = -t2642 * t2677 + t2643 * t2681;
t2556 = t2597 * t2638 + t2634 * t2601;
t2802 = t2547 * t2556;
t2777 = t2651 * t2665;
t2786 = t2650 * t2671;
t2579 = t2648 * t2777 - t2786;
t2752 = t2670 * t2649;
t2562 = t2664 * t2579 + t2648 * t2752;
t2801 = t2547 * t2562;
t2607 = t2634 * g(1) + t2638 * g(2);
t2800 = t2547 * t2607;
t2799 = t2547 * t2634;
t2798 = t2547 * t2638;
t2797 = t2648 * t2649;
t2795 = t2649 * t2653;
t2794 = t2649 * t2661;
t2793 = t2649 * t2663;
t2792 = t2649 * t2664;
t2791 = t2649 * t2665;
t2773 = t2652 * t2649;
t2772 = t2652 * t2653;
t2771 = t2652 * t2655;
t2770 = t2653 * t2654;
t2767 = t2660 * t2649;
t2766 = t2660 * t2661;
t2765 = t2660 * t2667;
t2764 = t2661 * t2666;
t2763 = t2662 * t2649;
t2762 = t2662 * t2663;
t2761 = t2662 * t2669;
t2760 = t2663 * t2668;
t2759 = t2664 * t2665;
t2758 = t2664 * t2671;
t2757 = t2665 * t2670;
t2742 = t2652 * t2825;
t2741 = t2654 * t2825;
t2740 = t2660 * t2824;
t2739 = t2666 * t2824;
t2738 = t2662 * t2823;
t2737 = t2668 * t2823;
t2736 = t2664 * t2822;
t2735 = t2670 * t2822;
t2573 = -t2648 * t2653 + t2650 * t2783;
t2533 = t2573 * t2829 + t2615 * t2790 - t2648 * t2846;
t2734 = t2533 * t2817;
t2733 = t2533 * t2814;
t2732 = t2533 * t2813;
t2580 = t2648 * t2661 - t2650 * t2776;
t2535 = t2580 * t2828 - t2619 * t2790 + t2648 * t2845;
t2731 = t2535 * t2809;
t2730 = t2535 * t2808;
t2729 = t2535 * t2812;
t2581 = t2648 * t2663 - t2650 * t2775;
t2537 = t2581 * t2827 - t2620 * t2790 + t2648 * t2844;
t2728 = t2537 * t2804;
t2727 = t2537 * t2803;
t2726 = t2537 * t2807;
t2582 = t2648 * t2665 - t2650 * t2774;
t2539 = t2582 * t2826 - t2621 * t2790 + t2648 * t2843;
t2725 = t2539 * t2799;
t2724 = t2539 * t2798;
t2723 = t2539 * t2802;
t2574 = t2648 * t2655 + t2650 * t2784;
t2561 = t2652 * t2574 + t2650 * t2769;
t2722 = t2561 * t2817;
t2721 = t2561 * t2814;
t2720 = t2561 * t2813;
t2583 = t2648 * t2667 + t2650 * t2781;
t2565 = t2660 * t2583 + t2650 * t2756;
t2719 = t2565 * t2812;
t2718 = t2565 * t2809;
t2717 = t2565 * t2808;
t2584 = t2648 * t2669 + t2650 * t2779;
t2566 = t2662 * t2584 + t2650 * t2754;
t2716 = t2566 * t2807;
t2715 = t2566 * t2804;
t2714 = t2566 * t2803;
t2585 = t2648 * t2671 + t2650 * t2777;
t2567 = t2585 * t2664 + t2650 * t2752;
t2713 = t2567 * t2802;
t2712 = t2567 * t2799;
t2711 = t2567 * t2798;
t2710 = t2561 * t2742;
t2709 = t2561 * t2741;
t2708 = t2565 * t2740;
t2707 = t2565 * t2739;
t2706 = t2566 * t2738;
t2705 = t2566 * t2737;
t2704 = t2567 * t2736;
t2703 = t2567 * t2735;
t2694 = t2746 * t2797;
t2693 = t2745 * t2797;
t2692 = t2744 * t2797;
t2691 = t2743 * t2797;
t2690 = pkin(3) * t2773 - t2651 * t2846;
t2689 = pkin(3) * t2767 - t2651 * t2845;
t2688 = pkin(3) * t2763 - t2651 * t2844;
t2687 = pkin(3) * t2792 - t2651 * t2843;
t2592 = t2650 * pkin(2) + t2672 * t2796;
t2593 = pkin(2) * t2796 - t2650 * t2672;
t2686 = t2572 * t2834 - (t2592 * t2655 - t2593 * t2653) * t2654;
t2685 = t2577 * t2833 - (t2592 * t2667 - t2593 * t2661) * t2666;
t2684 = t2578 * t2832 - (t2592 * t2669 - t2593 * t2663) * t2668;
t2683 = t2579 * t2831 - (t2592 * t2671 - t2593 * t2665) * t2670;
t2682 = 0.1e1 / pkin(3);
t2613 = t2643 * g(1) + t2642 * g(2);
t2612 = t2642 * g(1) - t2643 * g(2);
t2591 = t2651 * t2759 + t2752;
t2590 = t2651 * t2762 + t2754;
t2589 = t2651 * t2766 + t2756;
t2576 = t2651 * t2772 + t2769;
t2551 = t2621 * t2650 + t2687 * t2648;
t2550 = t2620 * t2650 + t2688 * t2648;
t2549 = t2619 * t2650 + t2689 * t2648;
t2548 = t2615 * t2650 + t2690 * t2648;
t2532 = -g(3) * t2585 - t2611 * t2579 + t2607 * t2791;
t2531 = -t2607 * t2649 * t2671 - g(3) * t2582 + t2611 * t2588;
t2530 = -t2606 * t2649 * t2669 - g(3) * t2581 + t2610 * t2587;
t2529 = -g(3) * t2584 - t2610 * t2578 + t2606 * t2793;
t2528 = -t2605 * t2649 * t2667 - g(3) * t2580 + t2609 * t2586;
t2527 = -g(3) * t2583 - t2609 * t2577 + t2605 * t2794;
t2526 = -t2604 * t2649 * t2655 + g(3) * t2573 + t2608 * t2575;
t2525 = -g(3) * t2574 - t2608 * t2572 + t2604 * t2795;
t2520 = (-t2589 * t2648 + t2650 * t2765) * t2609 - (t2589 * t2650 + t2648 * t2765) * g(3) - t2605 * (-t2649 * t2766 + t2651 * t2666);
t2519 = t2611 * (-t2591 * t2648 + t2650 * t2758) - g(3) * (t2591 * t2650 + t2648 * t2758) - t2607 * (-t2649 * t2759 + t2651 * t2670);
t2518 = t2610 * (-t2590 * t2648 + t2650 * t2761) - g(3) * (t2590 * t2650 + t2648 * t2761) - t2606 * (-t2649 * t2762 + t2651 * t2668);
t2517 = ((-t2651 * t2757 + t2792) * t2648 + t2670 * t2786) * t2611 + g(3) * (-t2585 * t2670 + t2650 * t2792) + t2607 * (t2665 * t2752 + t2778);
t2516 = ((-t2651 * t2760 + t2763) * t2648 + t2668 * t2787) * t2610 + g(3) * (-t2584 * t2668 + t2650 * t2763) + t2606 * (t2663 * t2754 + t2780);
t2515 = t2609 * ((-t2651 * t2764 + t2767) * t2648 + t2666 * t2788) + g(3) * (-t2583 * t2666 + t2650 * t2767) + t2605 * (t2661 * t2756 + t2782);
t2514 = (-t2576 * t2648 + t2650 * t2771) * t2608 - g(3) * (t2576 * t2650 + t2648 * t2771) - t2604 * (-t2649 * t2772 + t2651 * t2654);
t2513 = ((-t2651 * t2770 + t2773) * t2648 + t2654 * t2789) * t2608 + g(3) * (-t2574 * t2654 + t2650 * t2773) + t2604 * (t2653 * t2769 + t2785);
t1 = [-(-(t2579 * t2638 - t2634 * t2791) * t2831 + (t2551 * t2638 + t2571 * t2634) * t2670 + (t2651 * t2634 + t2638 * t2797) * t2835) * t2800 - (-(t2578 * t2637 - t2633 * t2793) * t2832 + (t2550 * t2637 + t2570 * t2633) * t2668 + (t2651 * t2633 + t2637 * t2797) * t2836) * t2805 - (-(t2577 * t2636 - t2632 * t2794) * t2833 + (t2549 * t2636 + t2569 * t2632) * t2666 + (t2651 * t2632 + t2636 * t2797) * t2837) * t2810 - ((-t2572 * t2635 + t2631 * t2795) * t2834 + (t2548 * t2635 + t2568 * t2631) * t2654 + (t2651 * t2631 + t2635 * t2797) * t2838) * t2815, 0, -t2526 * t2720 - t2528 * t2717 - t2530 * t2714 - t2531 * t2711, -t2525 * t2720 - t2527 * t2717 - t2529 * t2714 - t2532 * t2711, 0, 0, 0, 0, 0, -t2635 * t2709 - t2636 * t2707 - t2637 * t2705 - t2638 * t2703 + (-t2514 * t2732 + t2518 * t2727 + t2519 * t2724 + t2520 * t2730) * t2682, t2635 * t2710 + t2636 * t2708 + t2637 * t2706 + t2638 * t2704 + (-t2513 * t2732 + t2515 * t2730 + t2516 * t2727 + t2517 * t2724) * t2682, 0, 0, 0, -t2642 * t2612 - t2643 * t2613; -((t2579 * t2634 + t2638 * t2791) * t2831 + (-t2551 * t2634 + t2571 * t2638) * t2670 + (-t2634 * t2797 + t2651 * t2638) * t2835) * t2800 - ((t2578 * t2633 + t2637 * t2793) * t2832 + (-t2550 * t2633 + t2570 * t2637) * t2668 + (-t2633 * t2797 + t2651 * t2637) * t2836) * t2805 - ((t2577 * t2632 + t2636 * t2794) * t2833 + (-t2549 * t2632 + t2569 * t2636) * t2666 + (-t2632 * t2797 + t2651 * t2636) * t2837) * t2810 - (-(-t2572 * t2631 - t2635 * t2795) * t2834 + (-t2548 * t2631 + t2568 * t2635) * t2654 + (-t2631 * t2797 + t2651 * t2635) * t2838) * t2815, 0, t2526 * t2721 + t2528 * t2718 + t2530 * t2715 + t2531 * t2712, t2525 * t2721 + t2527 * t2718 + t2529 * t2715 + t2532 * t2712, 0, 0, 0, 0, 0, t2631 * t2709 + t2632 * t2707 + t2633 * t2705 + t2634 * t2703 + (t2514 * t2733 - t2518 * t2728 - t2519 * t2725 - t2520 * t2731) * t2682, -t2631 * t2710 - t2632 * t2708 - t2633 * t2706 - t2634 * t2704 + (t2513 * t2733 - t2515 * t2731 - t2516 * t2728 - t2517 * t2725) * t2682, 0, 0, 0, t2643 * t2612 - t2642 * t2613; -(-t2585 * t2831 - t2621 * t2648 * t2670 + (pkin(2) * t2792 + t2687 * t2670) * t2650) * t2800 - (-t2584 * t2832 - t2620 * t2648 * t2668 + (pkin(2) * t2763 + t2688 * t2668) * t2650) * t2805 - (-t2583 * t2833 - t2619 * t2648 * t2666 + (pkin(2) * t2767 + t2689 * t2666) * t2650) * t2810 - (-t2574 * t2834 - t2615 * t2648 * t2654 + (pkin(2) * t2773 + t2690 * t2654) * t2650) * t2815, 0, t2526 * t2816 + t2528 * t2811 + t2530 * t2806 + t2531 * t2801, t2525 * t2816 + t2527 * t2811 + t2529 * t2806 + t2532 * t2801, 0, 0, 0, 0, 0, t2560 * t2741 + t2563 * t2739 + t2564 * t2737 + t2562 * t2735 + (t2514 * t2821 + t2518 * t2819 + t2519 * t2818 + t2520 * t2820) * t2682, -t2560 * t2742 - t2563 * t2740 - t2564 * t2738 - t2562 * t2736 + (t2513 * t2821 + t2515 * t2820 + t2516 * t2819 + t2517 * t2818) * t2682, 0, 0, 0, -g(3); -((t2597 * t2683 + t2601 * t2850 + (-t2597 * t2691 + t2601 * t2700) * t2664) * t2638 + (-t2597 * t2850 + t2683 * t2601 + (-t2597 * t2700 - t2601 * t2691) * t2664) * t2634) / ((pkin(3) * t2757 + t2843) * t2752 + t2664 * t2700) * t2607 - ((t2596 * t2684 + t2600 * t2849 + (-t2596 * t2692 + t2600 * t2699) * t2662) * t2637 + (-t2596 * t2849 + t2684 * t2600 + (-t2596 * t2699 - t2600 * t2692) * t2662) * t2633) / ((pkin(3) * t2760 + t2844) * t2754 + t2662 * t2699) * t2606 - ((t2595 * t2685 + t2599 * t2848 + (-t2595 * t2693 + t2599 * t2701) * t2660) * t2636 + (-t2595 * t2848 + t2685 * t2599 + (-t2595 * t2701 - t2599 * t2693) * t2660) * t2632) / ((pkin(3) * t2764 + t2845) * t2756 + t2660 * t2701) * t2605 - ((t2686 * t2594 + t2598 * t2847 + (-t2594 * t2694 + t2598 * t2702) * t2652) * t2635 + (-t2594 * t2847 + t2598 * t2686 + (-t2594 * t2702 - t2598 * t2694) * t2652) * t2631) / ((pkin(3) * t2770 + t2846) * t2769 + t2652 * t2702) * t2604, 0, t2526 * t2722 + t2528 * t2719 + t2530 * t2716 + t2531 * t2713, t2525 * t2722 + t2527 * t2719 + t2529 * t2716 + t2532 * t2713, 0, 0, 0, 0, 0, t2553 * t2709 + t2554 * t2707 + t2555 * t2705 + t2556 * t2703 + (t2514 * t2734 - t2518 * t2726 - t2519 * t2723 - t2520 * t2729) * t2682, -t2553 * t2710 - t2554 * t2708 - t2555 * t2706 - t2556 * t2704 + (t2513 * t2734 - t2515 * t2729 - t2516 * t2726 - t2517 * t2723) * t2682, 0, t2612, t2613, 0;];
tau_reg  = t1;
