% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPRR1G3A0
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x8]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRR1G3A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G3A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3A0_coriolisvec_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:27:02
% EndTime: 2020-03-09 21:27:04
% DurationCPUTime: 2.74s
% Computational Cost: add. (34308->236), mult. (26523->369), div. (1473->10), fcn. (13608->59), ass. (0->202)
t794 = pkin(7) + qJ(3,3);
t802 = sin(qJ(3,3));
t720 = pkin(1) * sin(t794) + t802 * pkin(2);
t712 = 0.1e1 / t720 ^ 2;
t799 = legFrame(3,2);
t866 = qJ(1,3) + pkin(7);
t756 = t799 + t866;
t750 = qJ(3,3) + t756;
t757 = -t799 + t866;
t751 = qJ(3,3) + t757;
t693 = cos(t751) + cos(t750);
t762 = sin(qJ(1,3) + t794);
t815 = xDP(2);
t816 = xDP(1);
t872 = sin(t750) - sin(t751);
t814 = xDP(3);
t915 = -2 * t814;
t675 = t693 * t816 + t762 * t915 - t872 * t815;
t894 = t675 / 0.2e1;
t917 = t712 * t894;
t795 = pkin(7) + qJ(3,2);
t804 = sin(qJ(3,2));
t721 = pkin(1) * sin(t795) + t804 * pkin(2);
t715 = 0.1e1 / t721 ^ 2;
t800 = legFrame(2,2);
t867 = qJ(1,2) + pkin(7);
t758 = t800 + t867;
t752 = qJ(3,2) + t758;
t759 = -t800 + t867;
t753 = qJ(3,2) + t759;
t694 = cos(t753) + cos(t752);
t763 = sin(qJ(1,2) + t795);
t871 = sin(t752) - sin(t753);
t677 = t694 * t816 + t763 * t915 - t871 * t815;
t890 = t677 / 0.2e1;
t919 = t715 * t890;
t796 = pkin(7) + qJ(3,1);
t806 = sin(qJ(3,1));
t722 = pkin(1) * sin(t796) + t806 * pkin(2);
t718 = 0.1e1 / t722 ^ 2;
t801 = legFrame(1,2);
t868 = qJ(1,1) + pkin(7);
t760 = t801 + t868;
t754 = qJ(3,1) + t760;
t761 = -t801 + t868;
t755 = qJ(3,1) + t761;
t695 = cos(t755) + cos(t754);
t764 = sin(qJ(1,1) + t796);
t870 = sin(t754) - sin(t755);
t676 = t695 * t816 + t764 * t915 - t870 * t815;
t892 = t676 / 0.2e1;
t918 = t718 * t892;
t789 = qJ(1,2) + t800;
t790 = qJ(1,2) - t800;
t679 = t871 * pkin(3) + (sin(t758) - sin(t759)) * pkin(2) + (sin(t789) - sin(t790)) * pkin(1);
t682 = -t694 * pkin(3) + (-cos(t758) - cos(t759)) * pkin(2) + (-cos(t789) - cos(t790)) * pkin(1);
t688 = sin(qJ(1,2)) * pkin(1) + pkin(2) * sin(t867) + pkin(3) * t763;
t714 = 0.1e1 / t721;
t818 = 0.1e1 / pkin(3);
t864 = 2 * t814;
t899 = (t679 * t815 + t682 * t816 + t688 * t864) * t714 * t818;
t665 = t899 / 0.2e1;
t671 = t714 * t890;
t662 = t671 + t665;
t912 = pkin(3) * t662;
t650 = t714 * t665 * t912;
t785 = cos(t795);
t809 = cos(qJ(3,2));
t909 = t809 * pkin(2);
t651 = -t912 + (-pkin(1) * t785 - t909) * t671;
t656 = t671 + t665 / 0.2e1;
t798 = cos(pkin(7));
t797 = sin(pkin(7));
t876 = t797 * t804;
t684 = -t909 + (-t798 * t809 + t876) * pkin(1);
t777 = t798 * pkin(1) + pkin(2);
t914 = pkin(1) * t797;
t706 = t804 * t777 + t809 * t914;
t855 = t714 * t899;
t819 = pkin(1) ^ 2;
t793 = pkin(2) ^ 2 + t819;
t817 = pkin(3) ^ 2;
t869 = 0.2e1 * pkin(1);
t873 = 0.2e1 * pkin(2) * pkin(3);
t916 = pkin(2) * t798;
t902 = (t656 * t809 * t873 + t793 * t671 + t662 * t817 + (pkin(3) * t656 * t785 + t671 * t916) * t869) * t818;
t703 = -pkin(1) * t876 + t809 * t777;
t832 = (pkin(3) + t703) * t662 / t706 * t899;
t927 = -t832 / 0.2e1;
t940 = t656 * t706 * t855 + t714 * (0.2e1 * t650 + t927 + (-0.2e1 * t651 - t902) * t919) * t684;
t931 = t693 / 0.2e1;
t930 = t695 / 0.2e1;
t932 = -t870 / 0.2e1;
t933 = -t872 / 0.2e1;
t834 = (t665 + 0.2e1 * t671) * t684 * t855;
t882 = t706 * t714;
t845 = (t650 - t832 / 0.4e1 + (-t651 - t902 / 0.2e1) * t919) * t882;
t937 = t834 / 0.4e1 - t845;
t687 = sin(qJ(1,3)) * pkin(1) + pkin(2) * sin(t866) + pkin(3) * t762;
t936 = t687 / 0.4e1;
t935 = t688 / 0.4e1;
t689 = sin(qJ(1,1)) * pkin(1) + pkin(2) * sin(t868) + pkin(3) * t764;
t934 = t689 / 0.4e1;
t929 = t818 / 0.2e1;
t791 = qJ(1,1) + t801;
t792 = qJ(1,1) - t801;
t680 = t870 * pkin(3) + (sin(t760) - sin(t761)) * pkin(2) + (sin(t791) - sin(t792)) * pkin(1);
t683 = -t695 * pkin(3) + (-cos(t760) - cos(t761)) * pkin(2) + (-cos(t791) - cos(t792)) * pkin(1);
t717 = 0.1e1 / t722;
t900 = (t680 * t815 + t683 * t816 + t689 * t864) * t717 * t818;
t664 = t900 / 0.2e1;
t670 = t717 * t892;
t661 = t670 + t664;
t810 = cos(qJ(3,1));
t875 = t797 * t806;
t704 = -pkin(1) * t875 + t810 * t777;
t707 = t806 * t777 + t810 * t914;
t831 = (pkin(3) + t704) * t661 / t707 * t900;
t928 = -t831 / 0.2e1;
t787 = qJ(1,3) + t799;
t788 = qJ(1,3) - t799;
t678 = t872 * pkin(3) + (sin(t756) - sin(t757)) * pkin(2) + (sin(t787) - sin(t788)) * pkin(1);
t681 = -t693 * pkin(3) + (-cos(t756) - cos(t757)) * pkin(2) + (-cos(t787) - cos(t788)) * pkin(1);
t711 = 0.1e1 / t720;
t901 = (t678 * t815 + t681 * t816 + t687 * t864) * t711 * t818;
t663 = t901 / 0.2e1;
t669 = t711 * t894;
t660 = t669 + t663;
t808 = cos(qJ(3,3));
t877 = t797 * t802;
t702 = -pkin(1) * t877 + t808 * t777;
t705 = t802 * t777 + t808 * t914;
t833 = (pkin(3) + t702) * t660 / t705 * t901;
t926 = -t833 / 0.2e1;
t925 = t717 * t930;
t924 = t694 * t714 / 0.2e1;
t923 = t711 * t931;
t922 = t717 * t932;
t921 = -t871 * t714 / 0.2e1;
t920 = t711 * t933;
t913 = pkin(3) * t660;
t911 = t661 * pkin(3);
t910 = t808 * pkin(2);
t908 = t810 * pkin(2);
t636 = t650 + t927 + (-t651 - t902) * t919;
t907 = t636 * t714;
t648 = t717 * t664 * t911;
t786 = cos(t796);
t652 = -t911 + (-pkin(1) * t786 - t908) * t670;
t655 = t670 + t664 / 0.2e1;
t903 = (t655 * t810 * t873 + t793 * t670 + t661 * t817 + (pkin(3) * t655 * t786 + t670 * t916) * t869) * t818;
t638 = t648 + t928 + (-t652 - t903) * t918;
t906 = t638 * t717;
t649 = t711 * t663 * t913;
t784 = cos(t794);
t653 = -t913 + (-pkin(1) * t784 - t910) * t669;
t654 = t669 + t663 / 0.2e1;
t904 = (t654 * t808 * t873 + t793 * t669 + t660 * t817 + (pkin(3) * t654 * t784 + t669 * t916) * t869) * t818;
t640 = t649 + t926 + (-t653 - t904) * t917;
t905 = t640 * t711;
t898 = t675 ^ 2 * t711 * t712;
t897 = t676 ^ 2 * t717 * t718;
t896 = t677 ^ 2 * t714 * t715;
t883 = t705 * t711;
t881 = t707 * t717;
t880 = t711 * t762;
t879 = t714 * t763;
t878 = t717 * t764;
t874 = t818 / 0.8e1;
t645 = -t651 * t919 + t650;
t863 = t645 * t703 * t714;
t862 = t645 * t882;
t646 = -t652 * t918 + t648;
t861 = t646 * t704 * t717;
t860 = t646 * t881;
t647 = -t653 * t917 + t649;
t859 = t647 * t702 * t711;
t858 = t647 * t883;
t857 = t711 * t901;
t856 = t717 * t900;
t854 = t702 * t898;
t853 = t705 * t898;
t852 = t704 * t897;
t851 = t707 * t897;
t850 = t703 * t896;
t849 = t706 * t896;
t685 = t910 + (t798 * t808 - t877) * pkin(1);
t848 = t711 * (0.2e1 * t649 + t926 + (-0.2e1 * t653 - t904) * t917) * t685;
t847 = (t649 - t833 / 0.4e1 + (-t653 - t904 / 0.2e1) * t917) * t883;
t686 = t908 + (t798 * t810 - t875) * pkin(1);
t844 = t717 * (0.2e1 * t648 + t928 + (-0.2e1 * t652 - t903) * t918) * t686;
t843 = (t648 - t831 / 0.4e1 + (-t652 - t903 / 0.2e1) * t918) * t881;
t842 = -0.2e1 * t847;
t840 = -0.2e1 * t843;
t839 = t654 * t705 * t857;
t838 = t655 * t707 * t856;
t836 = (t663 + 0.2e1 * t669) * t685 * t857;
t835 = (t664 + 0.2e1 * t670) * t686 * t856;
t830 = -t645 * t879 - t646 * t878 - t647 * t880;
t829 = t645 * t921 + t646 * t922 + t647 * t920;
t828 = t645 * t924 + t646 * t925 + t647 * t923;
t1 = [t828, 0, 0, t819 * t828, t636 * t924 + t638 * t925 + t640 * t923 + (t681 * t905 + t682 * t907 + t683 * t906) * t929, (t681 * t853 + t682 * t849 + t683 * t851) * t874 - t693 * t839 / 0.2e1 - t695 * t838 / 0.2e1 + t844 * t930 + t848 * t931 + (t681 * t859 + t682 * t863 + t683 * t861) * t929 - t940 * t694 / 0.2e1, (t681 * t854 + t682 * t850 + t683 * t852) * t874 - t693 * t836 / 0.4e1 - t695 * t835 / 0.4e1 + t840 * t930 + t842 * t931 + (-t681 * t858 - t682 * t862 - t683 * t860) * t929 + t937 * t694, 0; t829, 0, 0, t819 * t829, t636 * t921 + t638 * t922 + t640 * t920 + (t678 * t905 + t679 * t907 + t680 * t906) * t929, (t678 * t853 + t679 * t849 + t680 * t851) * t874 + t838 * t870 / 0.2e1 + t839 * t872 / 0.2e1 + t844 * t932 + t848 * t933 + (t678 * t859 + t679 * t863 + t680 * t861) * t929 + t940 * t871 / 0.2e1, (t678 * t854 + t679 * t850 + t680 * t852) * t874 + t872 * t836 / 0.4e1 + t870 * t835 / 0.4e1 + t840 * t932 + t842 * t933 + (-t678 * t858 - t679 * t862 - t680 * t860) * t929 - t937 * t871, 0; t830, 0, 0, t830 * t819, -t636 * t879 - t638 * t878 - t640 * t880 + (t687 * t905 + t688 * t907 + t689 * t906) * t818, (t687 * t859 + t688 * t863 + t689 * t861 + t849 * t935 + t851 * t934 + t853 * t936) * t818 + (-t844 + t838) * t764 + t940 * t763 + (-t848 + t839) * t762, (-t687 * t858 - t688 * t862 - t689 * t860 + t850 * t935 + t852 * t934 + t854 * t936) * t818 + (0.2e1 * t843 + t835 / 0.2e1) * t764 + (0.2e1 * t845 - t834 / 0.2e1) * t763 + (0.2e1 * t847 + t836 / 0.2e1) * t762, 0;];
tau_reg  = t1;
