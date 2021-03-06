% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RRPRR8V2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [13x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR8V2G1A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 21:04
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR8V2G1A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(13,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_mdp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:03:57
% EndTime: 2020-08-06 21:04:10
% DurationCPUTime: 12.61s
% Computational Cost: add. (73596->475), mult. (127555->833), div. (12054->21), fcn. (70095->70), ass. (0->396)
t2639 = sin(qJ(2,3));
t2645 = cos(qJ(2,3));
t2630 = sin(pkin(7));
t2883 = pkin(3) * t2630;
t2631 = cos(pkin(7));
t2585 = t2631 * pkin(3);
t2924 = t2585 + pkin(2);
t2498 = t2639 * t2924 + t2645 * t2883;
t2654 = 0.2e1 * pkin(7);
t2663 = pkin(3) ^ 2;
t2890 = t2663 / 0.2e1;
t2517 = cos(t2654) * t2890 + pkin(2) * (t2585 + pkin(2) / 0.2e1);
t2652 = xDP(2);
t2602 = pkin(1) * t2652;
t2636 = pkin(5) + qJ(3,3);
t2604 = -pkin(6) - t2636;
t2653 = xDP(1);
t2528 = -t2604 * t2653 + t2602;
t2601 = t2653 * pkin(1);
t2531 = t2604 * t2652 + t2601;
t2633 = legFrame(3,3);
t2579 = t2633 + qJ(1,3);
t2656 = 0.2e1 * qJ(2,3);
t2617 = sin(t2656);
t2620 = cos(t2656);
t2825 = t2924 * t2645;
t2778 = 0.2e1 * t2825;
t2899 = 0.2e1 * t2517;
t2782 = t2653 * t2899;
t2783 = t2652 * t2899;
t2560 = pkin(2) * t2585;
t2665 = pkin(2) ^ 2;
t2603 = t2663 + t2665;
t2928 = 0.2e1 * t2560 + t2603;
t2797 = t2928 * t2653;
t2798 = t2928 * t2652;
t2819 = t2630 * t2639;
t2540 = pkin(2) * t2630 + sin(t2654) * pkin(3) / 0.2e1;
t2830 = t2540 * t2653;
t2831 = t2540 * t2652;
t2874 = 0.2e1 * pkin(3);
t2887 = pkin(3) * t2540;
t2651 = xDP(3);
t2926 = 2 * t2651;
t2459 = (t2620 * t2782 + t2531 * t2778 + (-t2531 * t2819 - t2617 * t2830) * t2874 + t2797) * cos(t2579) + (t2620 * t2783 + t2528 * t2778 + (-t2528 * t2819 - t2617 * t2831) * t2874 + t2798) * sin(t2579) + (pkin(1) * t2498 + t2517 * t2617 + t2620 * t2887) * t2926;
t2775 = pkin(3) * t2819;
t2495 = 0.1e1 / (-t2775 + t2825);
t2678 = t2604 ^ 2;
t2933 = t2495 * t2459 / t2678;
t2641 = sin(qJ(2,2));
t2647 = cos(qJ(2,2));
t2499 = t2641 * t2924 + t2647 * t2883;
t2637 = pkin(5) + qJ(3,2);
t2605 = -pkin(6) - t2637;
t2529 = -t2605 * t2653 + t2602;
t2532 = t2605 * t2652 + t2601;
t2634 = legFrame(2,3);
t2580 = t2634 + qJ(1,2);
t2658 = 0.2e1 * qJ(2,2);
t2618 = sin(t2658);
t2621 = cos(t2658);
t2824 = t2924 * t2647;
t2777 = 0.2e1 * t2824;
t2818 = t2630 * t2641;
t2460 = (t2621 * t2782 + t2532 * t2777 + (-t2532 * t2818 - t2618 * t2830) * t2874 + t2797) * cos(t2580) + (t2621 * t2783 + t2529 * t2777 + (-t2529 * t2818 - t2618 * t2831) * t2874 + t2798) * sin(t2580) + (pkin(1) * t2499 + t2517 * t2618 + t2621 * t2887) * t2926;
t2774 = pkin(3) * t2818;
t2496 = 0.1e1 / (-t2774 + t2824);
t2679 = t2605 ^ 2;
t2932 = t2496 * t2460 / t2679;
t2643 = sin(qJ(2,1));
t2649 = cos(qJ(2,1));
t2500 = t2643 * t2924 + t2649 * t2883;
t2638 = pkin(5) + qJ(3,1);
t2606 = -pkin(6) - t2638;
t2530 = -t2606 * t2653 + t2602;
t2533 = t2606 * t2652 + t2601;
t2635 = legFrame(1,3);
t2581 = t2635 + qJ(1,1);
t2660 = 0.2e1 * qJ(2,1);
t2619 = sin(t2660);
t2622 = cos(t2660);
t2823 = t2924 * t2649;
t2776 = 0.2e1 * t2823;
t2817 = t2630 * t2643;
t2461 = (t2622 * t2782 + t2533 * t2776 + (-t2533 * t2817 - t2619 * t2830) * t2874 + t2797) * cos(t2581) + (t2622 * t2783 + t2530 * t2776 + (-t2530 * t2817 - t2619 * t2831) * t2874 + t2798) * sin(t2581) + (pkin(1) * t2500 + t2517 * t2619 + t2622 * t2887) * t2926;
t2773 = pkin(3) * t2817;
t2497 = 0.1e1 / (-t2773 + t2823);
t2680 = t2606 ^ 2;
t2931 = t2497 * t2461 / t2680;
t2930 = 0.3e1 * t2665 + 0.6e1 * t2663;
t2929 = 0.6e1 * t2665 + 0.3e1 * t2663;
t2598 = t2645 * pkin(2);
t2612 = qJ(2,3) + pkin(7);
t2574 = cos(t2612);
t2886 = pkin(3) * t2574;
t2537 = t2598 + t2886;
t2599 = t2647 * pkin(2);
t2614 = qJ(2,2) + pkin(7);
t2576 = cos(t2614);
t2885 = pkin(3) * t2576;
t2538 = t2599 + t2885;
t2600 = t2649 * pkin(2);
t2616 = qJ(2,1) + pkin(7);
t2578 = cos(t2616);
t2884 = pkin(3) * t2578;
t2539 = t2600 + t2884;
t2927 = 0.2e1 * pkin(1);
t2925 = -pkin(5) / 0.2e1;
t2556 = t2598 + pkin(1);
t2557 = t2599 + pkin(1);
t2558 = t2600 + pkin(1);
t2586 = sin(t2633);
t2589 = cos(t2633);
t2640 = sin(qJ(1,3));
t2646 = cos(qJ(1,3));
t2510 = -t2586 * t2640 + t2589 * t2646;
t2513 = t2586 * t2646 + t2589 * t2640;
t2592 = 0.1e1 / t2604;
t2855 = t2495 * t2498;
t2480 = (t2510 * t2653 + t2513 * t2652 + t2651 * t2855) * t2592;
t2626 = t2645 ^ 2;
t2582 = 0.2e1 * t2626 - 0.1e1;
t2923 = t2480 * t2582;
t2587 = sin(t2634);
t2590 = cos(t2634);
t2642 = sin(qJ(1,2));
t2648 = cos(qJ(1,2));
t2511 = -t2587 * t2642 + t2590 * t2648;
t2514 = t2587 * t2648 + t2590 * t2642;
t2594 = 0.1e1 / t2605;
t2853 = t2496 * t2499;
t2481 = (t2511 * t2653 + t2514 * t2652 + t2651 * t2853) * t2594;
t2627 = t2647 ^ 2;
t2583 = 0.2e1 * t2627 - 0.1e1;
t2922 = t2481 * t2583;
t2588 = sin(t2635);
t2591 = cos(t2635);
t2644 = sin(qJ(1,1));
t2650 = cos(qJ(1,1));
t2512 = -t2588 * t2644 + t2591 * t2650;
t2515 = t2588 * t2650 + t2591 * t2644;
t2596 = 0.1e1 / t2606;
t2851 = t2497 * t2500;
t2482 = (t2512 * t2653 + t2515 * t2652 + t2651 * t2851) * t2596;
t2628 = t2649 ^ 2;
t2584 = 0.2e1 * t2628 - 0.1e1;
t2921 = t2482 * t2584;
t2917 = -0.8e1 * t2560 - 0.4e1 * t2603;
t2555 = 0.2e1 * t2616;
t2916 = t2665 * t2622 + t2663 * cos(t2555);
t2554 = 0.2e1 * t2614;
t2915 = t2665 * t2621 + t2663 * cos(t2554);
t2553 = 0.2e1 * t2612;
t2914 = t2665 * t2620 + t2663 * cos(t2553);
t2913 = t2619 * t2665 + sin(t2555) * t2663;
t2912 = t2618 * t2665 + sin(t2554) * t2663;
t2911 = t2617 * t2665 + sin(t2553) * t2663;
t2909 = 0.2e1 * pkin(5);
t2908 = 2 * MDP(5);
t2519 = 0.1e1 / t2537;
t2522 = 0.1e1 / t2538;
t2525 = 0.1e1 / t2539;
t2669 = t2537 ^ 2;
t2520 = 0.1e1 / t2669;
t2672 = t2538 ^ 2;
t2523 = 0.1e1 / t2672;
t2675 = t2539 ^ 2;
t2526 = 0.1e1 / t2675;
t2907 = 2 * MDP(12);
t2610 = t2631 ^ 2;
t2796 = t2560 + t2665 / 0.2e1;
t2900 = -0.2e1 * (t2610 - 0.1e1 / 0.2e1) * t2663 - 0.2e1 * t2796;
t2898 = -0.4e1 * pkin(1) * (t2890 + t2796);
t2664 = pkin(2) * t2665;
t2888 = pkin(2) * t2663;
t2897 = -0.2e1 * t2664 - 0.4e1 * t2888;
t2894 = pkin(2) * pkin(3);
t2893 = t2480 / 0.2e1;
t2892 = t2481 / 0.2e1;
t2891 = t2482 / 0.2e1;
t2629 = t2651 ^ 2;
t2889 = pkin(2) * t2629;
t2882 = pkin(3) * t2665;
t2881 = pkin(5) * t2629;
t2854 = t2495 * t2592;
t2766 = t2459 * t2854;
t2456 = t2766 / 0.2e1;
t2475 = pkin(1) * t2480;
t2450 = -t2475 + t2456;
t2763 = t2480 * t2854;
t2813 = t2663 * t2480;
t2829 = t2928 * t2629;
t2846 = t2520 * t2592;
t2426 = -(-t2480 * t2626 * t2900 + t2450 * t2775 - t2610 * t2813 + t2813 - (0.2e1 * t2480 * t2775 + t2450) * t2825) * t2763 - 0.1e1 / (t2598 + (t2631 * t2645 - t2819) * pkin(3)) * t2829 * t2846 + t2893 * t2933;
t2880 = t2426 * pkin(1);
t2852 = t2496 * t2594;
t2765 = t2460 * t2852;
t2457 = t2765 / 0.2e1;
t2476 = pkin(1) * t2481;
t2452 = -t2476 + t2457;
t2762 = t2481 * t2852;
t2812 = t2663 * t2481;
t2840 = t2523 * t2594;
t2427 = -(-t2481 * t2627 * t2900 + t2452 * t2774 - t2610 * t2812 + t2812 - (0.2e1 * t2481 * t2774 + t2452) * t2824) * t2762 - 0.1e1 / (t2599 + (t2631 * t2647 - t2818) * pkin(3)) * t2829 * t2840 + t2892 * t2932;
t2879 = t2427 * pkin(1);
t2850 = t2497 * t2596;
t2764 = t2461 * t2850;
t2458 = t2764 / 0.2e1;
t2474 = t2482 * pkin(1);
t2448 = -t2474 + t2458;
t2761 = t2482 * t2850;
t2811 = t2663 * t2482;
t2834 = t2526 * t2596;
t2428 = -(-t2482 * t2628 * t2900 + t2448 * t2773 - t2610 * t2811 + t2811 - (0.2e1 * t2482 * t2773 + t2448) * t2823) * t2761 - 0.1e1 / (t2600 + (t2631 * t2649 - t2817) * pkin(3)) * t2829 * t2834 + t2891 * t2931;
t2878 = t2428 * pkin(1);
t2877 = pkin(3) * t2927;
t2876 = -0.8e1 * t2894;
t2873 = t2426 * t2636;
t2872 = t2427 * t2637;
t2871 = t2428 * t2638;
t2847 = t2519 * t2651;
t2756 = t2639 * t2847;
t2727 = pkin(2) * t2756;
t2867 = (t2636 * t2893 + t2727) * t2480;
t2841 = t2522 * t2651;
t2751 = t2641 * t2841;
t2726 = pkin(2) * t2751;
t2866 = (t2637 * t2892 + t2726) * t2481;
t2835 = t2525 * t2651;
t2746 = t2643 * t2835;
t2725 = pkin(2) * t2746;
t2865 = (t2638 * t2891 + t2725) * t2482;
t2477 = t2480 ^ 2;
t2864 = t2477 * t2519;
t2478 = t2481 ^ 2;
t2863 = t2478 * t2522;
t2479 = t2482 ^ 2;
t2862 = t2479 * t2525;
t2861 = t2480 * t2639;
t2860 = t2480 * t2645;
t2859 = t2481 * t2641;
t2858 = t2481 * t2647;
t2857 = t2482 * t2643;
t2856 = t2482 * t2649;
t2849 = t2519 * t2639;
t2848 = t2519 * t2645;
t2845 = t2520 * t2645;
t2534 = pkin(2) * t2639 + pkin(3) * sin(t2612);
t2844 = t2519 * t2520 * t2534;
t2843 = t2522 * t2641;
t2842 = t2522 * t2647;
t2839 = t2523 * t2647;
t2535 = pkin(2) * t2641 + pkin(3) * sin(t2614);
t2838 = t2522 * t2523 * t2535;
t2837 = t2525 * t2643;
t2836 = t2525 * t2649;
t2833 = t2526 * t2649;
t2536 = pkin(2) * t2643 + pkin(3) * sin(t2616);
t2832 = t2525 * t2526 * t2536;
t2816 = t2639 * t2645;
t2815 = t2641 * t2647;
t2814 = t2643 * t2649;
t2754 = t2639 * t2844;
t2709 = t2629 * t2754;
t2706 = pkin(2) * t2709;
t2804 = t2459 * t2763 - t2845 * t2889 - t2706 + 0.2e1 * t2873;
t2749 = t2641 * t2838;
t2708 = t2629 * t2749;
t2705 = pkin(2) * t2708;
t2803 = t2460 * t2762 - t2839 * t2889 - t2705 + 0.2e1 * t2872;
t2744 = t2643 * t2832;
t2707 = t2629 * t2744;
t2704 = pkin(2) * t2707;
t2802 = t2461 * t2761 - t2833 * t2889 - t2704 + 0.2e1 * t2871;
t2666 = pkin(1) ^ 2;
t2795 = pkin(5) ^ 2 + t2666;
t2794 = t2645 * t2927;
t2793 = t2647 * t2927;
t2792 = t2649 * t2927;
t2791 = -0.3e1 * t2888;
t2790 = -0.2e1 * t2888;
t2789 = -0.3e1 * t2882;
t2788 = -0.2e1 * t2882;
t2611 = t2656 + pkin(7);
t2787 = sin(t2611) * t2894;
t2613 = t2658 + pkin(7);
t2786 = sin(t2613) * t2894;
t2615 = t2660 + pkin(7);
t2785 = sin(t2615) * t2894;
t2781 = -0.2e1 * t2847;
t2780 = -0.2e1 * t2841;
t2779 = -0.2e1 * t2835;
t2567 = cos(t2654 + qJ(2,3));
t2570 = cos(-pkin(7) + qJ(2,3));
t2573 = cos(t2611);
t2655 = 0.3e1 * qJ(2,3);
t2662 = pkin(3) * t2663;
t2769 = (-0.4e1 * t2537 * (pkin(1) * t2456 - (t2666 + t2678) * t2480) + (-0.4e1 * t2787 - 0.2e1 * t2911) * t2604 * t2847 + (t2662 * cos(0.3e1 * t2612) + t2886 * t2929 + t2664 * cos(t2655) + t2598 * t2930 - (cos(t2654 + t2655) + t2567) * t2791 - (cos(t2655 + pkin(7)) + t2570) * t2789) * t2480 + (t2573 * t2876 - 0.4e1 * t2914 + t2917) * (-t2475 + t2766 / 0.4e1)) * t2592 * t2480;
t2568 = cos(t2654 + qJ(2,2));
t2571 = cos(-pkin(7) + qJ(2,2));
t2575 = cos(t2613);
t2657 = 0.3e1 * qJ(2,2);
t2768 = (-0.4e1 * t2538 * (pkin(1) * t2457 - (t2666 + t2679) * t2481) + (-0.4e1 * t2786 - 0.2e1 * t2912) * t2605 * t2841 + (t2662 * cos(0.3e1 * t2614) + t2885 * t2929 + t2664 * cos(t2657) + t2599 * t2930 - (cos(t2654 + t2657) + t2568) * t2791 - (cos(t2657 + pkin(7)) + t2571) * t2789) * t2481 + (t2575 * t2876 - 0.4e1 * t2915 + t2917) * (-t2476 + t2765 / 0.4e1)) * t2594 * t2481;
t2569 = cos(t2654 + qJ(2,1));
t2572 = cos(-pkin(7) + qJ(2,1));
t2577 = cos(t2615);
t2659 = 0.3e1 * qJ(2,1);
t2767 = (-0.4e1 * t2539 * (pkin(1) * t2458 - (t2666 + t2680) * t2482) + (-0.4e1 * t2785 - 0.2e1 * t2913) * t2606 * t2835 + (t2662 * cos(0.3e1 * t2616) + t2884 * t2929 + t2664 * cos(t2659) + t2600 * t2930 - (cos(t2659 + t2654) + t2569) * t2791 - (cos(t2659 + pkin(7)) + t2572) * t2789) * t2482 + (t2577 * t2876 - 0.4e1 * t2916 + t2917) * (-t2474 + t2764 / 0.4e1)) * t2596 * t2482;
t2760 = t2498 * t2854;
t2759 = t2499 * t2852;
t2758 = t2500 * t2850;
t2757 = t2477 * t2849;
t2755 = t2645 * t2847;
t2753 = t2645 * t2844;
t2752 = t2478 * t2843;
t2750 = t2647 * t2841;
t2748 = t2647 * t2838;
t2747 = t2479 * t2837;
t2745 = t2649 * t2835;
t2743 = t2649 * t2832;
t2742 = t2847 / 0.2e1;
t2741 = t2841 / 0.2e1;
t2740 = t2835 / 0.2e1;
t2694 = 0.2e1 * t2787 + t2911;
t2728 = -0.2e1 * t2662 - 0.4e1 * t2882;
t2715 = (t2694 * t2604 * t2480 + (t2567 * t2790 + t2570 * t2788 + t2574 * t2728 + t2645 * t2897 + t2898) * t2847) * t2651 * t2846;
t2718 = (t2574 * t2877 + t2603 + (t2794 + (t2573 + t2631) * t2874) * pkin(2) + t2914) * t2480 * t2933;
t2737 = t2715 / 0.2e1 + (-t2769 / 0.4e1 + t2718 / 0.4e1) * t2519;
t2693 = 0.2e1 * t2786 + t2912;
t2714 = (t2693 * t2605 * t2481 + (t2568 * t2790 + t2571 * t2788 + t2576 * t2728 + t2647 * t2897 + t2898) * t2841) * t2651 * t2840;
t2717 = (t2576 * t2877 + t2603 + (t2793 + (t2575 + t2631) * t2874) * pkin(2) + t2915) * t2481 * t2932;
t2736 = t2714 / 0.2e1 + (-t2768 / 0.4e1 + t2717 / 0.4e1) * t2522;
t2692 = 0.2e1 * t2785 + t2913;
t2713 = (t2692 * t2606 * t2482 + (t2569 * t2790 + t2572 * t2788 + t2578 * t2728 + t2649 * t2897 + t2898) * t2835) * t2651 * t2834;
t2716 = (t2578 * t2877 + t2603 + (t2792 + (t2577 + t2631) * t2874) * pkin(2) + t2916) * t2482 * t2931;
t2735 = t2713 / 0.2e1 + (-t2767 / 0.4e1 + t2716 / 0.4e1) * t2525;
t2721 = t2426 * t2760;
t2720 = t2427 * t2759;
t2719 = t2428 * t2758;
t2712 = t2519 * t2760;
t2711 = t2522 * t2759;
t2710 = t2525 * t2758;
t2703 = t2760 * t2816;
t2702 = t2759 * t2815;
t2701 = t2758 * t2814;
t2700 = -t2520 * t2639 + t2753;
t2699 = t2754 + t2845;
t2698 = -t2523 * t2641 + t2748;
t2697 = t2749 + t2839;
t2696 = -t2526 * t2643 + t2743;
t2695 = t2744 + t2833;
t2408 = -t2556 * t2426 + t2737;
t2691 = -t2477 * MDP(11) + (t2408 - 0.2e1 * t2867) * MDP(12);
t2409 = -t2557 * t2427 + t2736;
t2690 = -t2478 * MDP(11) + (t2409 - 0.2e1 * t2866) * MDP(12);
t2410 = -t2558 * t2428 + t2735;
t2689 = -t2479 * MDP(11) + (t2410 - 0.2e1 * t2865) * MDP(12);
t2688 = t2426 * t2849 + t2427 * t2843 + t2428 * t2837;
t2687 = t2426 * t2848 + t2427 * t2842 + t2428 * t2836;
t2686 = t2592 * (t2699 * MDP(6) + t2700 * MDP(7));
t2685 = t2594 * (t2697 * MDP(6) + t2698 * MDP(7));
t2684 = t2596 * (t2695 * MDP(6) + t2696 * MDP(7));
t2405 = 0.2e1 * (t2880 - t2715 / 0.4e1 + (t2769 / 0.8e1 - t2718 / 0.8e1) * t2519) * t2598 - t2636 * t2706 - pkin(1) * t2737 + (t2626 * t2665 + (t2909 + qJ(3,3)) * qJ(3,3) + t2795) * t2426;
t2417 = -pkin(5) * t2709 + t2426 * t2794;
t2420 = -0.2e1 * t2639 * t2880 - t2753 * t2881;
t2435 = -pkin(2) * (-pkin(2) * t2861 + t2636 * t2742) * t2755 + t2480 * (pkin(1) * t2727 + t2456 * t2636);
t2462 = pkin(5) * t2645 * t2742 - pkin(1) * t2861;
t2465 = -pkin(1) * t2860 + t2756 * t2925;
t2623 = t2639 ^ 2;
t2683 = t2426 * MDP(1) + (t2426 * t2623 - 0.2e1 * t2755 * t2861) * MDP(4) + (t2462 * t2781 + t2417) * MDP(9) + (t2465 * t2781 + t2420) * MDP(10) + t2804 * MDP(11) + t2405 * MDP(12) + (t2426 * t2816 - t2847 * t2923) * t2908 + t2435 * t2907;
t2406 = 0.2e1 * (t2879 - t2714 / 0.4e1 + (t2768 / 0.8e1 - t2717 / 0.8e1) * t2522) * t2599 - t2637 * t2705 - pkin(1) * t2736 + (t2627 * t2665 + (t2909 + qJ(3,2)) * qJ(3,2) + t2795) * t2427;
t2418 = -pkin(5) * t2708 + t2427 * t2793;
t2421 = -0.2e1 * t2641 * t2879 - t2748 * t2881;
t2437 = -pkin(2) * (-pkin(2) * t2859 + t2637 * t2741) * t2750 + t2481 * (pkin(1) * t2726 + t2457 * t2637);
t2463 = pkin(5) * t2647 * t2741 - pkin(1) * t2859;
t2466 = -pkin(1) * t2858 + t2751 * t2925;
t2624 = t2641 ^ 2;
t2682 = t2427 * MDP(1) + (t2427 * t2624 - 0.2e1 * t2750 * t2859) * MDP(4) + (t2463 * t2780 + t2418) * MDP(9) + (t2466 * t2780 + t2421) * MDP(10) + t2803 * MDP(11) + t2406 * MDP(12) + (t2427 * t2815 - t2841 * t2922) * t2908 + t2437 * t2907;
t2407 = 0.2e1 * (t2878 - t2713 / 0.4e1 + (t2767 / 0.8e1 - t2716 / 0.8e1) * t2525) * t2600 - t2638 * t2704 - pkin(1) * t2735 + (t2628 * t2665 + (t2909 + qJ(3,1)) * qJ(3,1) + t2795) * t2428;
t2419 = -pkin(5) * t2707 + t2428 * t2792;
t2422 = -0.2e1 * t2643 * t2878 - t2743 * t2881;
t2436 = -pkin(2) * (-pkin(2) * t2857 + t2638 * t2740) * t2745 + t2482 * (pkin(1) * t2725 + t2458 * t2638);
t2464 = pkin(5) * t2649 * t2740 - pkin(1) * t2857;
t2467 = -pkin(1) * t2856 + t2746 * t2925;
t2625 = t2643 ^ 2;
t2681 = t2428 * MDP(1) + (t2428 * t2625 - 0.2e1 * t2745 * t2857) * MDP(4) + (t2464 * t2779 + t2419) * MDP(9) + (t2467 * t2779 + t2422) * MDP(10) + t2802 * MDP(11) + t2407 * MDP(12) + (t2428 * t2814 - t2835 * t2921) * t2908 + t2436 * t2907;
t2509 = t2557 * t2648 - t2605 * t2642;
t2508 = t2556 * t2646 - t2604 * t2640;
t2507 = t2558 * t2650 - t2606 * t2644;
t2506 = t2558 * t2644 + t2606 * t2650;
t2505 = t2557 * t2642 + t2605 * t2648;
t2504 = t2556 * t2640 + t2604 * t2646;
t2494 = t2536 * t2927 + t2692;
t2493 = t2535 * t2927 + t2693;
t2492 = t2534 * t2927 + t2694;
t1 = [-(t2689 * (-t2506 * t2588 + t2507 * t2591 + t2512 * t2884) + t2681 * t2512) * t2596 - (t2690 * (-t2505 * t2587 + t2509 * t2590 + t2511 * t2885) + t2682 * t2511) * t2594 - (t2691 * (-t2504 * t2586 + t2508 * t2589 + t2510 * t2886) + t2683 * t2510) * t2592 + (-t2510 * t2686 - t2511 * t2685 - t2512 * t2684) * t2629; -(t2689 * (t2506 * t2591 + t2507 * t2588 + t2515 * t2884) + t2681 * t2515) * t2596 - (t2690 * (t2505 * t2590 + t2509 * t2587 + t2514 * t2885) + t2682 * t2514) * t2594 - (t2691 * (t2504 * t2589 + t2508 * t2586 + t2513 * t2886) + t2683 * t2513) * t2592 + (-t2513 * t2686 - t2514 * t2685 - t2515 * t2684) * t2629; (-t2719 - t2720 - t2721) * MDP(1) + (-t2623 * t2721 - t2624 * t2720 - t2625 * t2719 - t2645 * t2757 - t2647 * t2752 - t2649 * t2747 + (t2480 * t2519 * t2703 + t2481 * t2522 * t2702 + t2482 * t2525 * t2701) * t2926) * MDP(4) + (-t2584 * t2862 - t2583 * t2863 - t2582 * t2864 - 0.2e1 * t2428 * t2701 - 0.2e1 * t2427 * t2702 - 0.2e1 * t2426 * t2703 + (t2710 * t2921 + t2711 * t2922 + t2712 * t2923) * t2926) * MDP(5) + ((-t2695 * t2758 - t2697 * t2759 - t2699 * t2760) * t2629 + t2688) * MDP(6) + ((-t2696 * t2758 - t2698 * t2759 - t2700 * t2760) * t2629 + t2687) * MDP(7) + (0.1e1 / t2675 ^ 2 * t2536 + 0.1e1 / t2672 ^ 2 * t2535 + 0.1e1 / t2669 ^ 2 * t2534) * t2629 * MDP(8) + (-t2417 * t2760 - t2418 * t2759 - t2419 * t2758 + (t2462 * t2712 + t2463 * t2711 + t2464 * t2710) * t2926 - t2688 * pkin(5) + (t2747 + t2752 + t2757) * pkin(1)) * MDP(9) + (-t2420 * t2760 - t2421 * t2759 - t2422 * t2758 + (t2465 * t2712 + t2466 * t2711 + t2467 * t2710) * t2926 - t2687 * pkin(5) + (t2477 * t2848 + t2478 * t2842 + t2479 * t2836) * pkin(1)) * MDP(10) + (-(-t2494 * t2862 / 0.2e1 + t2802 * t2851) * t2596 - (-t2493 * t2863 / 0.2e1 + t2803 * t2853) * t2594 - (-t2492 * t2864 / 0.2e1 + t2804 * t2855) * t2592 - t2688 * pkin(2)) * MDP(11) + (((t2832 * t2889 + (-t2871 - (-pkin(2) * t2856 - t2474 + t2764) * t2482) * t2643) * pkin(2) - (-t2865 + t2410 / 0.2e1) * t2596 * t2494) * t2525 + ((t2838 * t2889 + (-t2872 - (-pkin(2) * t2858 - t2476 + t2765) * t2481) * t2641) * pkin(2) - (-t2866 + t2409 / 0.2e1) * t2594 * t2493) * t2522 + ((t2844 * t2889 + (-t2873 - (-pkin(2) * t2860 - t2475 + t2766) * t2480) * t2639) * pkin(2) - (-t2867 + t2408 / 0.2e1) * t2592 * t2492) * t2519 + (-t2405 - 0.2e1 * t2435) * t2760 + (-t2406 - 0.2e1 * t2437) * t2759 + (-t2407 - 0.2e1 * t2436) * t2758) * MDP(12);];
taucX  = t1;
