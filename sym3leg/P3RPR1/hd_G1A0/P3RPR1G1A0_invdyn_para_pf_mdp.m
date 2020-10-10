% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [10x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPR1G1A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:58
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPR1G1A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(10,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1G1A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPR1G1A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPR1G1A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1G1A0_invdyn_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPR1G1A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1G1A0_invdyn_para_pf_mdp: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1G1A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1G1A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'P3RPR1G1A0_invdyn_para_pf_mdp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:58:14
% EndTime: 2019-05-03 14:58:16
% DurationCPUTime: 1.95s
% Computational Cost: add. (12837->266), mult. (20008->431), div. (648->9), fcn. (9609->14), ass. (0->182)
t800 = legFrame(3,3);
t788 = sin(t800);
t791 = cos(t800);
t756 = g(1) * t788 - g(2) * t791;
t759 = g(1) * t791 + g(2) * t788;
t806 = sin(qJ(1,3));
t809 = cos(qJ(1,3));
t719 = t756 * t809 + t759 * t806;
t801 = legFrame(2,3);
t789 = sin(t801);
t792 = cos(t801);
t757 = g(1) * t789 - g(2) * t792;
t760 = g(1) * t792 + g(2) * t789;
t807 = sin(qJ(1,2));
t810 = cos(qJ(1,2));
t720 = t757 * t810 + t760 * t807;
t802 = legFrame(1,3);
t790 = sin(t802);
t793 = cos(t802);
t758 = g(1) * t790 - g(2) * t793;
t761 = g(1) * t793 + g(2) * t790;
t808 = sin(qJ(1,1));
t811 = cos(qJ(1,1));
t721 = t758 * t811 + t761 * t808;
t887 = -t758 * t808 + t761 * t811;
t886 = -t757 * t807 + t760 * t810;
t885 = -t756 * t806 + t759 * t809;
t827 = 0.1e1 / qJ(2,1);
t823 = 0.1e1 / qJ(2,2);
t819 = 0.1e1 / qJ(2,3);
t884 = pkin(1) * g(2);
t815 = pkin(1) + pkin(2);
t830 = koppelP(3,2);
t833 = koppelP(3,1);
t774 = -qJ(2,3) * t830 + t815 * t833;
t777 = qJ(2,3) * t833 + t815 * t830;
t813 = xDP(2);
t786 = t813 * t815;
t817 = xP(3);
t796 = sin(t817);
t797 = cos(t817);
t812 = xDP(3);
t814 = xDP(1);
t701 = qJ(2,3) * t814 + t786 + (t774 * t797 - t777 * t796) * t812;
t787 = t814 * t815;
t704 = qJ(2,3) * t813 - t787 + (t774 * t796 + t777 * t797) * t812;
t686 = (t701 * t806 - t704 * t809) * t791 + (t701 * t809 + t704 * t806) * t788;
t750 = t796 * t833 + t797 * t830;
t731 = t750 * t812 - t814;
t753 = -t796 * t830 + t797 * t833;
t734 = t753 * t812 + t813;
t695 = (-t731 * t809 + t734 * t806) * t791 + t788 * (t731 * t806 + t734 * t809);
t798 = t812 ^ 2;
t803 = xDDP(3);
t804 = xDDP(2);
t707 = -t750 * t798 + t753 * t803 + t804;
t805 = xDDP(1);
t710 = -t750 * t803 - t753 * t798 + t805;
t739 = t788 * t809 + t791 * t806;
t740 = -t788 * t806 + t791 * t809;
t818 = qJ(2,3) ^ 2;
t821 = t819 / t818;
t877 = t695 * t686;
t856 = t821 * t877;
t820 = 0.1e1 / qJ(2,3) ^ 2;
t876 = t695 * t820;
t674 = -t856 + (t740 * t710 + t739 * t707 - (-t695 * t815 + t686) * t876) * t819;
t883 = pkin(1) * t674;
t831 = koppelP(2,2);
t834 = koppelP(2,1);
t775 = -qJ(2,2) * t831 + t815 * t834;
t778 = qJ(2,2) * t834 + t815 * t831;
t702 = qJ(2,2) * t814 + t786 + (t775 * t797 - t778 * t796) * t812;
t705 = qJ(2,2) * t813 - t787 + (t775 * t796 + t778 * t797) * t812;
t687 = (t702 * t807 - t705 * t810) * t792 + (t702 * t810 + t705 * t807) * t789;
t751 = t796 * t834 + t797 * t831;
t732 = t751 * t812 - t814;
t754 = -t796 * t831 + t797 * t834;
t735 = t754 * t812 + t813;
t696 = (-t732 * t810 + t735 * t807) * t792 + t789 * (t732 * t807 + t735 * t810);
t708 = -t751 * t798 + t754 * t803 + t804;
t711 = -t751 * t803 - t754 * t798 + t805;
t741 = t789 * t810 + t792 * t807;
t742 = -t789 * t807 + t792 * t810;
t822 = qJ(2,2) ^ 2;
t825 = t823 / t822;
t875 = t696 * t687;
t855 = t825 * t875;
t824 = 0.1e1 / qJ(2,2) ^ 2;
t874 = t696 * t824;
t675 = -t855 + (t742 * t711 + t741 * t708 - (-t696 * t815 + t687) * t874) * t823;
t882 = pkin(1) * t675;
t832 = koppelP(1,2);
t835 = koppelP(1,1);
t776 = -qJ(2,1) * t832 + t815 * t835;
t779 = qJ(2,1) * t835 + t815 * t832;
t703 = qJ(2,1) * t814 + t786 + (t776 * t797 - t779 * t796) * t812;
t706 = qJ(2,1) * t813 - t787 + (t776 * t796 + t779 * t797) * t812;
t688 = (t703 * t808 - t706 * t811) * t793 + (t703 * t811 + t706 * t808) * t790;
t752 = t796 * t835 + t797 * t832;
t733 = t752 * t812 - t814;
t755 = -t796 * t832 + t797 * t835;
t736 = t755 * t812 + t813;
t697 = (-t733 * t811 + t736 * t808) * t793 + t790 * (t733 * t808 + t736 * t811);
t709 = -t752 * t798 + t755 * t803 + t804;
t712 = -t752 * t803 - t755 * t798 + t805;
t743 = t790 * t811 + t793 * t808;
t744 = -t790 * t808 + t793 * t811;
t826 = qJ(2,1) ^ 2;
t829 = t827 / t826;
t873 = t697 * t688;
t854 = t829 * t873;
t828 = 0.1e1 / qJ(2,1) ^ 2;
t872 = t697 * t828;
t676 = -t854 + (t744 * t712 + t743 * t709 - (-t697 * t815 + t688) * t872) * t827;
t881 = pkin(1) * t676;
t692 = t695 ^ 2;
t880 = t692 * t821;
t693 = t696 ^ 2;
t879 = t693 * t825;
t694 = t697 ^ 2;
t878 = t694 * t829;
t853 = -pkin(2) ^ 2 + (-0.2e1 * pkin(2) - pkin(1)) * pkin(1);
t862 = ((-t818 + t853) * t695 + t686 * t815) * t819 * t876 + t815 * t856;
t861 = ((-t822 + t853) * t696 + t687 * t815) * t823 * t874 + t815 * t855;
t860 = ((-t826 + t853) * t697 + t688 * t815) * t827 * t872 + t815 * t854;
t859 = 0.2e1 * t877;
t858 = 0.2e1 * t875;
t857 = 0.2e1 * t873;
t852 = t862 + t883;
t851 = t861 + t882;
t850 = t860 + t881;
t768 = -qJ(2,3) * t809 + t806 * t815;
t771 = qJ(2,3) * t806 + t809 * t815;
t713 = t768 * t791 + t771 * t788;
t716 = -t768 * t788 + t771 * t791;
t846 = -t707 * t713 - t710 * t716;
t849 = -MDP(4) * t674 + MDP(6) * ((-t692 - t846) * t819 - t852 - t719);
t769 = -qJ(2,2) * t810 + t807 * t815;
t772 = qJ(2,2) * t807 + t810 * t815;
t714 = t769 * t792 + t772 * t789;
t717 = -t769 * t789 + t772 * t792;
t845 = -t708 * t714 - t711 * t717;
t848 = -MDP(4) * t675 + MDP(6) * ((-t693 - t845) * t823 - t851 - t720);
t770 = -qJ(2,1) * t811 + t808 * t815;
t773 = qJ(2,1) * t808 + t811 * t815;
t715 = t770 * t793 + t773 * t790;
t718 = -t770 * t790 + t773 * t793;
t844 = -t709 * t715 - t712 * t718;
t847 = -MDP(4) * t676 + MDP(6) * ((-t694 - t844) * t827 - t850 - t721);
t816 = pkin(1) * g(1);
t780 = -g(2) * qJ(2,3) + t816;
t781 = g(1) * qJ(2,3) + t884;
t843 = MDP(1) * t674 + MDP(2) * t719 + MDP(3) * t885 + MDP(4) * (t819 * t846 + t719 + t862 + 0.2e1 * t883) + MDP(5) * (0.2e1 * qJ(2,3) * t674 + t820 * t859 - t885) + MDP(6) * ((t780 * t806 - t781 * t809) * t791 + (t780 * t809 + t781 * t806) * t788 + t674 * t818 + t852 * pkin(1) + (pkin(1) * t846 + t859) * t819);
t782 = -g(2) * qJ(2,2) + t816;
t783 = g(1) * qJ(2,2) + t884;
t842 = MDP(1) * t675 + MDP(2) * t720 + MDP(3) * t886 + MDP(4) * (t823 * t845 + t720 + t861 + 0.2e1 * t882) + MDP(5) * (0.2e1 * qJ(2,2) * t675 + t824 * t858 - t886) + MDP(6) * ((t782 * t807 - t783 * t810) * t792 + (t782 * t810 + t783 * t807) * t789 + t675 * t822 + t851 * pkin(1) + (pkin(1) * t845 + t858) * t823);
t784 = -g(2) * qJ(2,1) + t816;
t785 = g(1) * qJ(2,1) + t884;
t841 = MDP(1) * t676 + MDP(2) * t721 + MDP(3) * t887 + MDP(4) * (t827 * t844 + t721 + t860 + 0.2e1 * t881) + MDP(5) * (0.2e1 * qJ(2,1) * t676 + t828 * t857 - t887) + MDP(6) * ((t784 * t808 - t785 * t811) * t793 + (t784 * t811 + t785 * t808) * t790 + t676 * t826 + t850 * pkin(1) + (pkin(1) * t844 + t857) * t827);
t795 = t805 - g(1);
t794 = t804 - g(2);
t767 = t808 * t832 + t811 * t835;
t766 = t808 * t835 - t811 * t832;
t765 = t807 * t831 + t810 * t834;
t764 = t807 * t834 - t810 * t831;
t763 = t806 * t830 + t809 * t833;
t762 = t806 * t833 - t809 * t830;
t749 = -t796 * t803 - t797 * t798;
t748 = -t796 * t798 + t797 * t803;
t738 = t794 * t796 + t795 * t797;
t737 = t794 * t797 - t795 * t796;
t730 = t776 * t808 - t779 * t811;
t729 = t775 * t807 - t778 * t810;
t728 = t774 * t806 - t777 * t809;
t727 = t776 * t811 + t779 * t808;
t726 = t775 * t810 + t778 * t807;
t725 = t774 * t809 + t777 * t806;
t691 = (-t727 * t796 + t730 * t797) * t793 + (t727 * t797 + t730 * t796) * t790;
t690 = (-t726 * t796 + t729 * t797) * t792 + (t726 * t797 + t729 * t796) * t789;
t689 = (-t725 * t796 + t728 * t797) * t791 + (t725 * t797 + t728 * t796) * t788;
t1 = [(-t716 * t880 - t717 * t879 - t718 * t878) * MDP(5) + t749 * MDP(8) - t748 * MDP(9) + (-t737 * t796 + t738 * t797) * MDP(10) + (t718 * t847 + t744 * t841) * t827 + (t717 * t848 + t742 * t842) * t823 + (t716 * t849 + t740 * t843) * t819; (-t713 * t880 - t714 * t879 - t715 * t878) * MDP(5) + t748 * MDP(8) + t749 * MDP(9) + (t737 * t797 + t738 * t796) * MDP(10) + (t715 * t847 + t743 * t841) * t827 + (t714 * t848 + t741 * t842) * t823 + (t713 * t849 + t739 * t843) * t819; (-t689 * t880 - t690 * t879 - t691 * t878) * MDP(5) + t803 * MDP(7) + t737 * MDP(8) - t738 * MDP(9) + (t847 * t691 + t841 * ((t766 * t797 - t767 * t796) * t793 + t790 * (t766 * t796 + t767 * t797))) * t827 + (t848 * t690 + t842 * ((t764 * t797 - t765 * t796) * t792 + t789 * (t764 * t796 + t765 * t797))) * t823 + (t849 * t689 + t843 * ((t762 * t797 - t763 * t796) * t791 + t788 * (t762 * t796 + t763 * t797))) * t819;];
tauX  = t1;
