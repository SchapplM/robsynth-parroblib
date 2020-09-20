% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRR1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:23
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRR1G1P1A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:20
% EndTime: 2020-03-09 21:23:21
% DurationCPUTime: 1.53s
% Computational Cost: add. (17915->249), mult. (12308->402), div. (1908->11), fcn. (7980->44), ass. (0->193)
t827 = 1 / pkin(3);
t884 = 2 * t827;
t883 = MDP(4) * pkin(1);
t802 = pkin(7) + qJ(3,3);
t785 = sin(t802);
t882 = pkin(1) * t785;
t803 = pkin(7) + qJ(3,2);
t786 = sin(t803);
t881 = pkin(1) * t786;
t804 = pkin(7) + qJ(3,1);
t787 = sin(t804);
t880 = pkin(1) * t787;
t788 = cos(t802);
t879 = pkin(1) * t788;
t789 = cos(t803);
t878 = pkin(1) * t789;
t790 = cos(t804);
t877 = pkin(1) * t790;
t876 = pkin(1) * sin(pkin(7));
t806 = cos(pkin(7));
t875 = pkin(2) * t806;
t809 = legFrame(1,3);
t793 = t809 + qJ(1,1);
t777 = pkin(7) + t793;
t774 = qJ(3,1) + t777;
t762 = sin(t774);
t765 = cos(t774);
t824 = xDP(2);
t825 = xDP(1);
t723 = t762 * t824 + t765 * t825;
t816 = sin(qJ(3,1));
t750 = pkin(2) * t816 + t880;
t746 = 0.1e1 / t750;
t705 = t723 * t746;
t714 = -pkin(1) * sin(t793) - pkin(2) * sin(t777) - pkin(3) * t762;
t717 = -pkin(1) * cos(t793) - pkin(2) * cos(t777) - pkin(3) * t765;
t696 = t714 * t824 + t717 * t825;
t841 = t696 * t827 * t746;
t690 = t705 + t841 / 0.2e1;
t693 = t705 + t841;
t801 = pkin(1) ^ 2 + pkin(2) ^ 2;
t822 = cos(qJ(3,1));
t826 = pkin(3) ^ 2;
t844 = 0.2e1 * pkin(1);
t848 = 0.2e1 * pkin(2) * pkin(3);
t874 = (t690 * t822 * t848 + t801 * t705 + t693 * t826 + (pkin(3) * t690 * t790 + t705 * t875) * t844) * t827;
t807 = legFrame(3,3);
t791 = t807 + qJ(1,3);
t775 = pkin(7) + t791;
t772 = qJ(3,3) + t775;
t760 = sin(t772);
t763 = cos(t772);
t721 = t760 * t824 + t763 * t825;
t812 = sin(qJ(3,3));
t748 = pkin(2) * t812 + t882;
t742 = 0.1e1 / t748;
t703 = t721 * t742;
t712 = -pkin(1) * sin(t791) - pkin(2) * sin(t775) - pkin(3) * t760;
t715 = -pkin(1) * cos(t791) - pkin(2) * cos(t775) - pkin(3) * t763;
t695 = t712 * t824 + t715 * t825;
t842 = t695 * t827 * t742;
t688 = t703 + t842 / 0.2e1;
t692 = t703 + t842;
t818 = cos(qJ(3,3));
t873 = (t688 * t818 * t848 + t801 * t703 + t692 * t826 + (pkin(3) * t688 * t788 + t703 * t875) * t844) * t827;
t808 = legFrame(2,3);
t792 = t808 + qJ(1,2);
t776 = pkin(7) + t792;
t773 = qJ(3,2) + t776;
t761 = sin(t773);
t764 = cos(t773);
t722 = t761 * t824 + t764 * t825;
t814 = sin(qJ(3,2));
t749 = pkin(2) * t814 + t881;
t744 = 0.1e1 / t749;
t704 = t722 * t744;
t713 = -pkin(1) * sin(t792) - pkin(2) * sin(t776) - pkin(3) * t761;
t716 = -pkin(1) * cos(t792) - pkin(2) * cos(t776) - pkin(3) * t764;
t694 = t713 * t824 + t716 * t825;
t843 = t694 * t827 * t744;
t689 = t704 + t843 / 0.2e1;
t691 = t704 + t843;
t820 = cos(qJ(3,2));
t872 = (t689 * t820 * t848 + t801 * t704 + t691 * t826 + (pkin(3) * t689 * t789 + t704 * t875) * t844) * t827;
t871 = t712 * t742;
t870 = t713 * t744;
t869 = t714 * t746;
t868 = t715 * t742;
t867 = t716 * t744;
t866 = t717 * t746;
t718 = t721 ^ 2;
t743 = 0.1e1 / t748 ^ 2;
t865 = t718 * t743;
t719 = t722 ^ 2;
t745 = 0.1e1 / t749 ^ 2;
t864 = t719 * t745;
t720 = t723 ^ 2;
t747 = 0.1e1 / t750 ^ 2;
t863 = t720 * t747;
t862 = t721 * t743;
t861 = t722 * t745;
t860 = t723 * t747;
t859 = t742 * t760;
t858 = t742 * t763;
t857 = t743 * t818;
t856 = t744 * t761;
t855 = t744 * t764;
t854 = t745 * t820;
t853 = t746 * t762;
t852 = t746 * t765;
t851 = t747 * t822;
t810 = xDDP(2);
t850 = t810 * t827;
t811 = xDDP(1);
t849 = t811 * t827;
t847 = g(1) * t763 + g(2) * t760;
t846 = g(1) * t764 + g(2) * t761;
t845 = g(1) * t765 + g(2) * t762;
t685 = t691 * t745 * t694;
t725 = t810 * t856;
t728 = t811 * t855;
t840 = t685 + t725 + t728;
t686 = t692 * t743 * t695;
t724 = t810 * t859;
t727 = t811 * t858;
t839 = t686 + t724 + t727;
t687 = t693 * t747 * t696;
t726 = t810 * t853;
t729 = t811 * t852;
t838 = t687 + t726 + t729;
t837 = g(1) * t760 - g(2) * t763;
t836 = g(1) * t761 - g(2) * t764;
t835 = g(1) * t762 - g(2) * t765;
t834 = -0.2e1 * t688 * t842;
t833 = -0.2e1 * t689 * t843;
t832 = -0.2e1 * t690 * t841;
t784 = pkin(1) * t806 + pkin(2);
t831 = -t691 * (t784 * t820 - t814 * t876 + pkin(3)) / (t784 * t814 + t820 * t876) * t843 + t849 * t867 + t850 * t870;
t830 = -t692 * (t784 * t818 - t812 * t876 + pkin(3)) / (t784 * t812 + t818 * t876) * t842 + t849 * t868 + t850 * t871;
t829 = -t693 * (t784 * t822 - t816 * t876 + pkin(3)) / (t784 * t816 + t822 * t876) * t841 + t849 * t866 + t850 * t869;
t828 = 0.1e1 / pkin(3) ^ 2;
t823 = cos(qJ(1,1));
t821 = cos(qJ(1,2));
t819 = cos(qJ(1,3));
t817 = sin(qJ(1,1));
t815 = sin(qJ(1,2));
t813 = sin(qJ(1,3));
t799 = cos(t809);
t798 = cos(t808);
t797 = cos(t807);
t796 = sin(t809);
t795 = sin(t808);
t794 = sin(t807);
t741 = g(1) * t799 + g(2) * t796;
t740 = g(1) * t798 + g(2) * t795;
t739 = g(1) * t797 + g(2) * t794;
t738 = g(1) * t796 - g(2) * t799;
t737 = g(1) * t795 - g(2) * t798;
t736 = g(1) * t794 - g(2) * t797;
t711 = -t738 * t817 + t741 * t823;
t710 = -t737 * t815 + t740 * t821;
t709 = -t736 * t813 + t739 * t819;
t708 = t738 * t823 + t741 * t817;
t707 = t737 * t821 + t740 * t815;
t706 = t736 * t819 + t739 * t813;
t684 = -pkin(3) * t691 + (-pkin(2) * t820 - t878) * t704;
t683 = -t693 * pkin(3) + (-pkin(2) * t822 - t877) * t705;
t682 = -pkin(3) * t692 + (-pkin(2) * t818 - t879) * t703;
t678 = -t683 * t860 + t838;
t677 = -t684 * t861 + t840;
t676 = -t682 * t862 + t839;
t675 = (g(1) * t817 - g(2) * t823) * t799 + (g(1) * t823 + g(2) * t817) * t796 + pkin(1) * t678;
t674 = (g(1) * t815 - g(2) * t821) * t798 + (g(1) * t821 + g(2) * t815) * t795 + pkin(1) * t677;
t673 = (g(1) * t813 - g(2) * t819) * t797 + (g(1) * t819 + g(2) * t813) * t794 + pkin(1) * t676;
t672 = (-t678 * t816 + t720 * t851) * pkin(2) + (-t678 * t787 + t790 * t863) * pkin(1) + t845;
t671 = (-t677 * t814 + t719 * t854) * pkin(2) + (-t677 * t786 + t789 * t864) * pkin(1) + t846;
t670 = (t676 * t818 + t812 * t865) * pkin(2) + (t676 * t788 + t785 * t865) * pkin(1) + t837;
t669 = (t677 * t820 + t814 * t864) * pkin(2) + (t677 * t789 + t786 * t864) * pkin(1) + t836;
t668 = (-t676 * t812 + t718 * t857) * pkin(2) + (-t676 * t785 + t788 * t865) * pkin(1) + t847;
t667 = (t678 * t822 + t816 * t863) * pkin(2) + (t678 * t790 + t787 * t863) * pkin(1) + t835;
t666 = (-t683 - t874) * t860 + t829 + t838;
t665 = (-t684 - t872) * t861 + t831 + t840;
t664 = (-t682 - t873) * t862 + t830 + t839;
t663 = 0.2e1 * t687 + 0.2e1 * t726 + 0.2e1 * t729 + (-0.2e1 * t683 - t874) * t860 + t829;
t662 = 0.2e1 * t685 + 0.2e1 * t725 + 0.2e1 * t728 + (-0.2e1 * t684 - t872) * t861 + t831;
t661 = 0.2e1 * t686 + 0.2e1 * t724 + 0.2e1 * t727 + (-0.2e1 * t682 - t873) * t862 + t830;
t660 = t832 * t877 - t663 * t880 - pkin(2) * (t816 * t663 + (t696 * t828 + t723 * t884) * t696 * t851) + t845;
t659 = t833 * t878 - t662 * t881 - pkin(2) * (t814 * t662 + (t694 * t828 + t722 * t884) * t694 * t854) + t846;
t658 = t834 * t879 - t661 * t882 - pkin(2) * (t812 * t661 + (t695 * t828 + t721 * t884) * t695 * t857) + t847;
t657 = pkin(2) * (t663 * t822 + t816 * t832) + (t663 * t790 + t787 * t832) * pkin(1) + t835;
t656 = pkin(2) * (t662 * t820 + t814 * t833) + (t662 * t789 + t786 * t833) * pkin(1) + t836;
t655 = pkin(2) * (t661 * t818 + t812 * t834) + (t661 * t788 + t785 * t834) * pkin(1) + t837;
t1 = [(t676 * t858 + t677 * t855 + t678 * t852) * MDP(1) + (t706 * t858 + t707 * t855 + t708 * t852) * MDP(2) + (t709 * t858 + t710 * t855 + t711 * t852) * MDP(3) + (t664 * t858 + t665 * t855 + t666 * t852) * MDP(5) + (t655 * t858 + t656 * t855 + t657 * t852) * MDP(6) + (t658 * t858 + t659 * t855 + t660 * t852) * MDP(7) + (t811 - g(1)) * MDP(8) + (t673 * t858 + t674 * t855 + t675 * t852) * t883 + ((t664 * t868 + t665 * t867 + t666 * t866) * MDP(5) + (t667 * t866 + t669 * t867 + t670 * t868) * MDP(6) + (t668 * t868 + t671 * t867 + t672 * t866) * MDP(7)) * t827; (t676 * t859 + t677 * t856 + t678 * t853) * MDP(1) + (t706 * t859 + t707 * t856 + t708 * t853) * MDP(2) + (t709 * t859 + t710 * t856 + t711 * t853) * MDP(3) + (t664 * t859 + t665 * t856 + t666 * t853) * MDP(5) + (t655 * t859 + t656 * t856 + t657 * t853) * MDP(6) + (t658 * t859 + t659 * t856 + t660 * t853) * MDP(7) + (t810 - g(2)) * MDP(8) + (t673 * t859 + t674 * t856 + t675 * t853) * t883 + ((t664 * t871 + t665 * t870 + t666 * t869) * MDP(5) + (t667 * t869 + t669 * t870 + t670 * t871) * MDP(6) + (t668 * t871 + t671 * t870 + t672 * t869) * MDP(7)) * t827; ((3 * MDP(4)) + MDP(8)) * (xDDP(3) - g(3));];
tauX  = t1;
