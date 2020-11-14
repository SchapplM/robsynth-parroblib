% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR8V2G1A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR8V2G1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V2G1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_mdp: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:36:32
% EndTime: 2020-08-06 17:36:34
% DurationCPUTime: 1.63s
% Computational Cost: add. (871->208), mult. (1839->405), div. (54->7), fcn. (2022->22), ass. (0->161)
t2748 = cos(qJ(2,1));
t2742 = sin(qJ(2,1));
t2749 = pkin(7) + pkin(6);
t2769 = t2742 * t2749;
t2712 = pkin(2) * t2748 + t2769;
t2730 = sin(pkin(8));
t2732 = cos(pkin(8));
t2720 = t2749 * t2748;
t2709 = pkin(2) * t2742 - t2720;
t2733 = cos(pkin(4));
t2731 = sin(pkin(4));
t2741 = sin(qJ(3,1));
t2772 = t2741 * t2731;
t2751 = pkin(3) * t2772 - t2709 * t2733;
t2824 = t2712 * t2732 + t2751 * t2730;
t2746 = cos(qJ(2,2));
t2740 = sin(qJ(2,2));
t2773 = t2740 * t2749;
t2711 = pkin(2) * t2746 + t2773;
t2719 = t2749 * t2746;
t2708 = pkin(2) * t2740 - t2719;
t2739 = sin(qJ(3,2));
t2776 = t2739 * t2731;
t2752 = pkin(3) * t2776 - t2708 * t2733;
t2823 = t2711 * t2732 + t2752 * t2730;
t2744 = cos(qJ(2,3));
t2738 = sin(qJ(2,3));
t2777 = t2738 * t2749;
t2710 = pkin(2) * t2744 + t2777;
t2718 = t2749 * t2744;
t2707 = pkin(2) * t2738 - t2718;
t2737 = sin(qJ(3,3));
t2780 = t2737 * t2731;
t2753 = pkin(3) * t2780 - t2707 * t2733;
t2822 = t2710 * t2732 + t2753 * t2730;
t2821 = MDP(1) * g(3);
t2743 = cos(qJ(3,3));
t2820 = pkin(3) * t2743 ^ 2;
t2745 = cos(qJ(3,2));
t2819 = pkin(3) * t2745 ^ 2;
t2747 = cos(qJ(3,1));
t2818 = pkin(3) * t2747 ^ 2;
t2817 = t2731 * g(3);
t2792 = t2733 * t2737;
t2795 = t2731 * t2738;
t2668 = 0.1e1 / (t2795 * t2820 + (pkin(3) * t2792 + t2707 * t2731) * t2743 + pkin(2) * t2792);
t2713 = t2730 * g(1) - t2732 * g(2);
t2714 = t2732 * g(1) + t2730 * g(2);
t2734 = legFrame(3,3);
t2721 = sin(t2734);
t2724 = cos(t2734);
t2816 = ((-t2817 + (t2713 * t2724 + t2714 * t2721) * t2733) * t2744 + t2738 * (-t2721 * t2713 + t2714 * t2724)) * t2668;
t2790 = t2733 * t2739;
t2794 = t2731 * t2740;
t2669 = 0.1e1 / (t2794 * t2819 + (pkin(3) * t2790 + t2708 * t2731) * t2745 + pkin(2) * t2790);
t2735 = legFrame(2,3);
t2722 = sin(t2735);
t2725 = cos(t2735);
t2815 = ((-t2817 + (t2713 * t2725 + t2714 * t2722) * t2733) * t2746 + t2740 * (-t2722 * t2713 + t2714 * t2725)) * t2669;
t2788 = t2733 * t2741;
t2793 = t2731 * t2742;
t2670 = 0.1e1 / (t2793 * t2818 + (pkin(3) * t2788 + t2709 * t2731) * t2747 + pkin(2) * t2788);
t2736 = legFrame(1,3);
t2723 = sin(t2736);
t2726 = cos(t2736);
t2814 = ((-t2817 + (t2713 * t2726 + t2714 * t2723) * t2733) * t2748 + (-t2723 * t2713 + t2714 * t2726) * t2742) * t2670;
t2715 = t2743 * pkin(3) + pkin(2);
t2698 = t2738 * t2715 - t2718;
t2768 = t2743 * t2731;
t2674 = 0.1e1 / (t2698 * t2768 + t2715 * t2792);
t2677 = -t2730 * t2721 + t2724 * t2732;
t2680 = t2732 * t2721 + t2724 * t2730;
t2801 = (t2715 * t2744 + t2777) * t2733;
t2813 = (-t2698 * t2677 - t2680 * t2801) * t2674;
t2716 = t2745 * pkin(3) + pkin(2);
t2699 = t2740 * t2716 - t2719;
t2766 = t2745 * t2731;
t2675 = 0.1e1 / (t2699 * t2766 + t2716 * t2790);
t2678 = -t2730 * t2722 + t2725 * t2732;
t2681 = t2732 * t2722 + t2725 * t2730;
t2800 = (t2716 * t2746 + t2773) * t2733;
t2812 = (-t2699 * t2678 - t2681 * t2800) * t2675;
t2717 = t2747 * pkin(3) + pkin(2);
t2700 = t2742 * t2717 - t2720;
t2764 = t2747 * t2731;
t2676 = 0.1e1 / (t2700 * t2764 + t2717 * t2788);
t2679 = -t2730 * t2723 + t2726 * t2732;
t2682 = t2732 * t2723 + t2726 * t2730;
t2799 = (t2717 * t2748 + t2769) * t2733;
t2811 = (-t2700 * t2679 - t2682 * t2799) * t2676;
t2810 = (-t2677 * t2801 + t2698 * t2680) * t2674;
t2809 = (-t2678 * t2800 + t2699 * t2681) * t2675;
t2808 = (-t2679 * t2799 + t2700 * t2682) * t2676;
t2701 = t2721 * g(1) - t2724 * g(2);
t2704 = t2724 * g(1) + t2721 * g(2);
t2785 = t2733 * t2744;
t2807 = (t2704 * (t2730 * t2785 + t2732 * t2738) + t2701 * (-t2730 * t2738 + t2732 * t2785) - t2744 * t2817) * t2668;
t2791 = t2733 * t2738;
t2683 = t2730 * t2791 - t2732 * t2744;
t2686 = t2730 * t2744 + t2732 * t2791;
t2806 = (g(3) * t2795 - t2704 * t2683 - t2701 * t2686) * t2668;
t2702 = t2722 * g(1) - t2725 * g(2);
t2705 = t2725 * g(1) + t2722 * g(2);
t2783 = t2733 * t2746;
t2805 = (t2705 * (t2730 * t2783 + t2732 * t2740) + t2702 * (-t2730 * t2740 + t2732 * t2783) - t2746 * t2817) * t2669;
t2789 = t2733 * t2740;
t2684 = t2730 * t2789 - t2732 * t2746;
t2687 = t2730 * t2746 + t2732 * t2789;
t2804 = (g(3) * t2794 - t2705 * t2684 - t2702 * t2687) * t2669;
t2703 = t2723 * g(1) - t2726 * g(2);
t2706 = t2726 * g(1) + t2723 * g(2);
t2781 = t2733 * t2748;
t2803 = (t2706 * (t2730 * t2781 + t2732 * t2742) + t2703 * (-t2730 * t2742 + t2732 * t2781) - t2748 * t2817) * t2670;
t2787 = t2733 * t2742;
t2685 = t2730 * t2787 - t2732 * t2748;
t2688 = t2730 * t2748 + t2732 * t2787;
t2802 = (g(3) * t2793 - t2706 * t2685 - t2703 * t2688) * t2670;
t2786 = t2733 * t2743;
t2784 = t2733 * t2745;
t2782 = t2733 * t2747;
t2779 = t2737 * t2738;
t2778 = t2737 * t2744;
t2775 = t2739 * t2740;
t2774 = t2739 * t2746;
t2771 = t2741 * t2742;
t2770 = t2741 * t2748;
t2767 = t2743 * t2744;
t2765 = t2745 * t2746;
t2763 = t2747 * t2748;
t2762 = pkin(2) * t2780;
t2761 = pkin(2) * t2776;
t2760 = pkin(2) * t2772;
t2759 = t2737 * t2816;
t2758 = t2743 * t2816;
t2757 = t2739 * t2815;
t2756 = t2745 * t2815;
t2755 = t2741 * t2814;
t2754 = t2747 * t2814;
t2750 = 0.1e1 / pkin(3);
t2694 = t2742 * t2782 - t2772;
t2693 = t2740 * t2784 - t2776;
t2692 = t2738 * t2786 - t2780;
t2691 = t2733 * t2771 + t2764;
t2690 = t2733 * t2775 + t2766;
t2689 = t2733 * t2779 + t2768;
t2673 = t2730 * t2712 - t2751 * t2732;
t2672 = t2730 * t2711 - t2752 * t2732;
t2671 = t2730 * t2710 - t2753 * t2732;
t2655 = -t2682 * t2764 - (-t2748 * t2679 + t2682 * t2787) * t2741;
t2654 = -t2681 * t2766 - (-t2746 * t2678 + t2681 * t2789) * t2739;
t2653 = -t2680 * t2768 - (-t2744 * t2677 + t2680 * t2791) * t2737;
t2652 = -t2679 * t2764 - (t2679 * t2787 + t2748 * t2682) * t2741;
t2651 = -t2678 * t2766 - (t2678 * t2789 + t2746 * t2681) * t2739;
t2650 = -t2677 * t2768 - (t2677 * t2791 + t2744 * t2680) * t2737;
t2646 = (-t2694 * t2730 + t2732 * t2763) * t2706 - t2703 * (t2694 * t2732 + t2730 * t2763) + g(3) * (t2742 * t2764 + t2788);
t2645 = (-t2693 * t2730 + t2732 * t2765) * t2705 - t2702 * (t2693 * t2732 + t2730 * t2765) + g(3) * (t2740 * t2766 + t2790);
t2644 = (-t2692 * t2730 + t2732 * t2767) * t2704 - t2701 * (t2692 * t2732 + t2730 * t2767) + g(3) * (t2738 * t2768 + t2792);
t2643 = t2706 * (-t2691 * t2730 + t2732 * t2770) - (t2691 * t2732 + t2730 * t2770) * t2703 - g(3) * (-t2731 * t2771 + t2782);
t2642 = t2705 * (-t2690 * t2730 + t2732 * t2774) - (t2690 * t2732 + t2730 * t2774) * t2702 - g(3) * (-t2731 * t2775 + t2784);
t2641 = t2704 * (-t2689 * t2730 + t2732 * t2778) - (t2689 * t2732 + t2730 * t2778) * t2701 - g(3) * (-t2731 * t2779 + t2786);
t1 = [(-(-(t2685 * t2726 + t2723 * t2688) * t2818 + (-t2673 * t2723 + t2824 * t2726) * t2747 + t2682 * t2760) * t2670 - (-(t2684 * t2725 + t2722 * t2687) * t2819 + (-t2672 * t2722 + t2823 * t2725) * t2745 + t2681 * t2761) * t2669 - (-(t2683 * t2724 + t2721 * t2686) * t2820 + (-t2671 * t2721 + t2822 * t2724) * t2743 + t2680 * t2762) * t2668) * t2821 + (t2650 * t2807 + t2651 * t2805 + t2652 * t2803) * MDP(3) + (t2650 * t2806 + t2651 * t2804 + t2652 * t2802) * MDP(4) + (t2650 * t2758 + t2651 * t2756 + t2652 * t2754 + (t2641 * t2810 + t2642 * t2809 + t2643 * t2808) * t2750) * MDP(10) + (-t2650 * t2759 - t2651 * t2757 - t2652 * t2755 + (t2644 * t2810 + t2645 * t2809 + t2646 * t2808) * t2750) * MDP(11) - g(1) * MDP(12); (-((-t2723 * t2685 + t2688 * t2726) * t2818 + (t2673 * t2726 + t2824 * t2723) * t2747 - t2679 * t2760) * t2670 - ((-t2722 * t2684 + t2687 * t2725) * t2819 + (t2672 * t2725 + t2823 * t2722) * t2745 - t2678 * t2761) * t2669 - ((-t2721 * t2683 + t2686 * t2724) * t2820 + (t2671 * t2724 + t2822 * t2721) * t2743 - t2677 * t2762) * t2668) * t2821 + (t2653 * t2807 + t2654 * t2805 + t2655 * t2803) * MDP(3) + (t2653 * t2806 + t2654 * t2804 + t2655 * t2802) * MDP(4) + (t2653 * t2758 + t2654 * t2756 + t2655 * t2754 + (t2641 * t2813 + t2642 * t2812 + t2643 * t2811) * t2750) * MDP(10) + (-t2653 * t2759 - t2654 * t2757 - t2655 * t2755 + (t2644 * t2813 + t2645 * t2812 + t2646 * t2811) * t2750) * MDP(11) - g(2) * MDP(12); (-(3 * MDP(1)) - MDP(12)) * g(3);];
taugX  = t1;
