% Calculate vector of centrifugal and Coriolis load for parallel robot
% P3RRPRR8V2G2A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-07 13:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:06:02
% EndTime: 2022-11-07 13:06:07
% DurationCPUTime: 5.90s
% Computational Cost: add. (26973->439), mult. (46998->672), div. (4179->12), fcn. (33657->50), ass. (0->316)
t815 = qJ(3,1) + pkin(5);
t748 = t815 * mrSges(3,2) - Ifges(3,6);
t810 = sin(pkin(7));
t811 = cos(pkin(7));
t993 = t815 * mrSges(3,1) - Ifges(3,5);
t1009 = -t748 * t810 + t993 * t811;
t814 = qJ(3,2) + pkin(5);
t747 = t814 * mrSges(3,2) - Ifges(3,6);
t994 = t814 * mrSges(3,1) - Ifges(3,5);
t1008 = -t747 * t810 + t994 * t811;
t813 = qJ(3,3) + pkin(5);
t746 = t813 * mrSges(3,2) - Ifges(3,6);
t995 = t813 * mrSges(3,1) - Ifges(3,5);
t1007 = -t746 * t810 + t995 * t811;
t816 = Ifges(3,2) - Ifges(3,1);
t740 = -pkin(2) * mrSges(3,2) - t816 * t810;
t949 = t810 * mrSges(3,1);
t745 = -pkin(2) * t949 + Ifges(2,4) - Ifges(3,4);
t798 = t811 ^ 2;
t951 = Ifges(3,4) * t798;
t1006 = -t740 * t811 - t745 - 0.2e1 * t951;
t1005 = -0.2e1 * pkin(2);
t910 = 2 * pkin(3);
t1004 = pkin(1) / 0.2e1;
t1003 = pkin(2) / 0.2e1;
t804 = qJ(2,1) + pkin(7);
t773 = cos(t804);
t758 = pkin(3) * t773;
t832 = cos(qJ(2,1));
t790 = t832 * pkin(2);
t912 = t790 + t758;
t1002 = pkin(1) + t912;
t803 = qJ(2,2) + pkin(7);
t772 = cos(t803);
t757 = pkin(3) * t772;
t830 = cos(qJ(2,2));
t789 = t830 * pkin(2);
t913 = t789 + t757;
t1001 = pkin(1) + t913;
t802 = qJ(2,3) + pkin(7);
t771 = cos(t802);
t756 = pkin(3) * t771;
t828 = cos(qJ(2,3));
t788 = t828 * pkin(2);
t914 = t788 + t756;
t1000 = pkin(1) + t914;
t992 = 0.2e1 * pkin(1);
t991 = 0.2e1 * mrSges(3,1);
t990 = -2 * mrSges(3,3);
t732 = 0.1e1 / t914;
t819 = legFrame(3,2);
t778 = sin(t819);
t781 = cos(t819);
t837 = xDP(2);
t838 = xDP(1);
t702 = (t778 * t838 + t781 * t837) * t732;
t699 = t702 ^ 2;
t733 = 0.1e1 / t913;
t820 = legFrame(2,2);
t779 = sin(t820);
t782 = cos(t820);
t703 = (t779 * t838 + t782 * t837) * t733;
t700 = t703 ^ 2;
t734 = 0.1e1 / t912;
t821 = legFrame(1,2);
t780 = sin(t821);
t783 = cos(t821);
t704 = (t780 * t838 + t783 * t837) * t734;
t701 = t704 ^ 2;
t956 = t810 * pkin(3);
t767 = pkin(1) * t956;
t822 = sin(qJ(2,3));
t856 = pkin(3) ^ 2;
t857 = pkin(2) ^ 2;
t955 = t811 * pkin(3);
t766 = pkin(2) * t955;
t981 = 0.2e1 * t766;
t864 = 0.2e1 * t798 * t856 - t856 + t857 + t981;
t715 = t864 * t822 + t767;
t823 = sin(qJ(1,3));
t795 = pkin(6) + t813;
t829 = cos(qJ(1,3));
t882 = pkin(1) * t823 - t829 * t795;
t927 = t810 * t822;
t897 = pkin(3) * t927;
t722 = -0.2e1 * t823 * t897 + t882;
t812 = t857 / 0.2e1;
t915 = t766 + t812;
t729 = (t798 - 0.1e1 / 0.2e1) * t856 + t915;
t963 = pkin(1) * t822;
t741 = -t956 + t963;
t759 = pkin(2) + t955;
t903 = t781 * t956;
t935 = t778 * t823;
t936 = t778 * t759;
t939 = t759 * t781;
t894 = pkin(3) * (t798 - 0.1e1);
t959 = (t823 * t894 + t882 * t927) * pkin(3);
t806 = t828 ^ 2;
t980 = 0.2e1 * t806;
t678 = (-t729 * t935 + t759 * t903) * t980 + (t781 * t715 - t722 * t936) * t828 + t778 * t959 + t741 * t939;
t874 = -t759 * t828 + t897;
t725 = 0.1e1 / t874;
t775 = 0.1e1 / t795;
t943 = t725 * t775;
t669 = t678 * t837 * t943;
t900 = t778 * t956;
t930 = t781 * t823;
t679 = (t729 * t930 + t759 * t900) * t980 + (t715 * t778 + t722 * t939) * t828 - t781 * t959 + t741 * t936;
t670 = t679 * t838 * t943;
t719 = t1000 * t829 + t823 * t795;
t836 = xDP(3);
t709 = t719 * t775 * t836;
t663 = -t670 - t669 + t709;
t989 = 0.2e1 * t663;
t824 = sin(qJ(2,2));
t716 = t864 * t824 + t767;
t825 = sin(qJ(1,2));
t796 = pkin(6) + t814;
t831 = cos(qJ(1,2));
t881 = pkin(1) * t825 - t831 * t796;
t926 = t810 * t824;
t896 = pkin(3) * t926;
t723 = -0.2e1 * t825 * t896 + t881;
t962 = pkin(1) * t824;
t742 = -t956 + t962;
t902 = t782 * t956;
t933 = t779 * t825;
t934 = t779 * t759;
t938 = t759 * t782;
t958 = (t825 * t894 + t881 * t926) * pkin(3);
t807 = t830 ^ 2;
t979 = 0.2e1 * t807;
t680 = (-t729 * t933 + t759 * t902) * t979 + (t782 * t716 - t723 * t934) * t830 + t779 * t958 + t742 * t938;
t873 = -t759 * t830 + t896;
t726 = 0.1e1 / t873;
t776 = 0.1e1 / t796;
t942 = t726 * t776;
t671 = t680 * t837 * t942;
t899 = t779 * t956;
t929 = t782 * t825;
t681 = (t729 * t929 + t759 * t899) * t979 + (t716 * t779 + t723 * t938) * t830 - t782 * t958 + t742 * t934;
t672 = t681 * t838 * t942;
t720 = t1001 * t831 + t825 * t796;
t710 = t720 * t776 * t836;
t664 = -t672 - t671 + t710;
t988 = 0.2e1 * t664;
t826 = sin(qJ(2,1));
t717 = t864 * t826 + t767;
t827 = sin(qJ(1,1));
t797 = pkin(6) + t815;
t833 = cos(qJ(1,1));
t880 = pkin(1) * t827 - t833 * t797;
t925 = t810 * t826;
t895 = pkin(3) * t925;
t724 = -0.2e1 * t827 * t895 + t880;
t961 = pkin(1) * t826;
t743 = -t956 + t961;
t901 = t783 * t956;
t931 = t780 * t827;
t932 = t780 * t759;
t937 = t759 * t783;
t957 = (t827 * t894 + t880 * t925) * pkin(3);
t808 = t832 ^ 2;
t978 = 0.2e1 * t808;
t682 = (-t729 * t931 + t759 * t901) * t978 + (t783 * t717 - t724 * t932) * t832 + t780 * t957 + t743 * t937;
t872 = -t759 * t832 + t895;
t727 = 0.1e1 / t872;
t777 = 0.1e1 / t797;
t941 = t727 * t777;
t673 = t682 * t837 * t941;
t898 = t780 * t956;
t928 = t783 * t827;
t683 = (t729 * t928 + t759 * t898) * t978 + (t717 * t780 + t724 * t937) * t832 - t783 * t957 + t743 * t932;
t674 = t683 * t838 * t941;
t721 = t1002 * t833 + t827 * t797;
t711 = t721 * t777 * t836;
t665 = -t674 - t673 + t711;
t987 = 0.2e1 * t665;
t693 = (-t759 * t935 + t903) * t828 + (t823 * t900 + t939) * t822;
t696 = (t759 * t930 + t900) * t828 + t822 * (-t823 * t903 + t936);
t687 = (t829 * t836 - (t693 * t837 + t696 * t838) * t725) * t775;
t986 = 0.2e1 * t687;
t694 = (-t759 * t933 + t902) * t830 + (t825 * t899 + t938) * t824;
t697 = (t759 * t929 + t899) * t830 + t824 * (-t825 * t902 + t934);
t688 = (t831 * t836 - (t694 * t837 + t697 * t838) * t726) * t776;
t985 = 0.2e1 * t688;
t695 = (-t759 * t931 + t901) * t832 + (t827 * t898 + t937) * t826;
t698 = (t759 * t928 + t898) * t832 + t826 * (-t827 * t901 + t932);
t689 = (t833 * t836 - (t695 * t837 + t698 * t838) * t727) * t777;
t984 = 0.2e1 * t689;
t924 = t816 * t798;
t950 = Ifges(3,4) * t810;
t952 = mrSges(3,2) * t810;
t965 = m(3) * t857;
t708 = 0.2e1 * t924 + (pkin(2) * t991 + 0.4e1 * t950) * t811 + t965 + t952 * t1005 + Ifges(2,2) - Ifges(2,1) - t816;
t983 = -0.2e1 * t708;
t982 = -0.4e1 * pkin(1) * (t856 / 0.2e1 + t915);
t977 = -0.2e1 * t822;
t976 = -0.2e1 * t824;
t975 = -0.2e1 * t826;
t974 = -0.4e1 * t828;
t973 = -0.4e1 * t830;
t972 = -0.4e1 * t832;
t971 = -4 * pkin(5) - 4 * pkin(6);
t844 = m(3) * pkin(2);
t970 = mrSges(2,2) * pkin(1);
t969 = pkin(5) * mrSges(2,2);
t968 = t732 / 0.2e1;
t967 = t733 / 0.2e1;
t966 = t734 / 0.2e1;
t784 = mrSges(2,1) + t844;
t964 = pkin(1) * t784;
t960 = pkin(3) * t857;
t684 = t689 * pkin(1);
t954 = (mrSges(2,3) + mrSges(3,3)) * pkin(5);
t953 = t856 * pkin(2);
t948 = pkin(5) * mrSges(2,1) - Ifges(2,5);
t947 = -t784 * pkin(5) + Ifges(2,5);
t946 = t1006 * t806;
t945 = t1006 * t807;
t944 = t1006 * t808;
t923 = t822 * t828;
t922 = t824 * t830;
t921 = t826 * t832;
t920 = -0.2e1 * t844;
t919 = pkin(2) * t910;
t685 = pkin(1) * t687;
t658 = t685 + t670 / 0.2e1 + t669 / 0.2e1 - t709 / 0.2e1;
t736 = pkin(2) * t822 + pkin(3) * sin(t802);
t760 = 0.2e1 * t802;
t765 = pkin(2) * t812 + t953;
t851 = 0.2e1 * qJ(2,3);
t801 = t851 + pkin(7);
t849 = 0.2e1 * pkin(7);
t854 = pkin(5) ^ 2;
t858 = pkin(1) ^ 2;
t911 = t856 + t857;
t871 = -(2 * t854) - 0.2e1 * t858 - t911 + ((-4 * pkin(5) - 2 * pkin(6)) * pkin(6));
t886 = -(2 * pkin(3) * t856) - 0.4e1 * t960;
t907 = -0.2e1 * t953;
t909 = -0.2e1 * t960;
t642 = ((cos(qJ(2,3) - pkin(7)) * t909 + cos(t849 + qJ(2,3)) * t907 + t886 * t771 + t765 * t974 + t982) * t699 * t968 + (t1000 * t663 - 0.2e1 * t658 * t756 + t989 * t1004 + (-t856 * cos(t760) - t857 * cos(t851) + t871 + ((t971 - 2 * qJ(3,3)) * qJ(3,3))) * t687 / 0.2e1 + (t658 * t974 + (-cos(t801) - t811) * t687 * t910) * t1003 + (t736 + (sin(t801) * t919 + sin(t760) * t856 + sin(t851) * t857) * t968) * t702 * t795) * t687) * t775;
t744 = t981 + t911;
t654 = (-t699 * t744 / (t788 + (t811 * t828 - t927) * pkin(3)) + (t874 * t687 - t685 + t989) * t687) * t775;
t675 = t687 * t964;
t887 = -m(3) * qJ(3,3) - mrSges(3,3);
t890 = Ifges(2,6) - t969;
t690 = -(t887 * pkin(2) - t1007 + t947) * t822 - (-t746 * t811 - t810 * t995 + t890) * t828;
t879 = mrSges(3,1) * t811 - t952;
t730 = -t844 - t879;
t845 = m(3) * pkin(1);
t878 = mrSges(3,2) * t811 + t949;
t705 = -t730 * t828 - t878 * t822 + t845;
t735 = 0.2e1 * t740;
t739 = (t784 - t952) * t992;
t843 = m(3) * pkin(5);
t752 = t843 - t887;
t763 = -t969 / 0.2e1 + Ifges(2,6) / 0.2e1;
t850 = -pkin(5) / 0.2e1;
t785 = -qJ(3,3) / 0.2e1 + t850;
t805 = pkin(1) * t991;
t841 = Ifges(3,6) / 0.2e1;
t842 = Ifges(3,5) / 0.2e1;
t861 = -(m(2) * t854) - (m(2) + m(3)) * t858 - Ifges(2,1) - Ifges(3,2) - Ifges(1,3) + t924;
t862 = 0.2e1 * t1006;
t863 = 0.2e1 * mrSges(2,2) + 0.2e1 * t878;
t870 = t879 * t822;
t893 = t699 * t732 * t736;
t906 = 0.4e1 * t951;
t908 = (mrSges(2,2) + t949) * t992;
t918 = (-t708 * t806 - (t822 * t906 + (t735 * t822 + t805) * t811 + t739) * t828 + t822 * t908 - t813 ^ 2 * m(3) + (qJ(3,3) * t990) + t861) * t654 - t690 * t893 + t705 * t642 + 0.2e1 * (-t745 * t923 - (-mrSges(3,2) * t963 - t950) * t811 - t954) * t654 + t663 * t752 * t986 + (t675 * t977 + (-(t752 * pkin(2) + t1007 + t948) * t828 + ((t785 * mrSges(3,2) + t841) * t811 + (t785 * mrSges(3,1) + t842) * t810 + t763) * t977) * t702 + (t923 * t983 - 0.4e1 * t946 + (-t828 * t863 - 0.2e1 * t870) * pkin(1) + t862) * t687) * t702;
t686 = pkin(1) * t688;
t659 = t686 + t672 / 0.2e1 + t671 / 0.2e1 - t710 / 0.2e1;
t738 = pkin(2) * t824 + pkin(3) * sin(t803);
t761 = 0.2e1 * t803;
t852 = 0.2e1 * qJ(2,2);
t799 = pkin(7) + t852;
t643 = ((cos(qJ(2,2) - pkin(7)) * t909 + cos(t849 + qJ(2,2)) * t907 + t886 * t772 + t765 * t973 + t982) * t700 * t967 + (t1001 * t664 - 0.2e1 * t659 * t757 + t988 * t1004 + (-t857 * cos(t852) - t856 * cos(t761) + t871 + ((t971 - 2 * qJ(3,2)) * qJ(3,2))) * t688 / 0.2e1 + (t659 * t973 + (-cos(t799) - t811) * t688 * t910) * t1003 + (t738 + (sin(t799) * t919 + sin(t761) * t856 + sin(t852) * t857) * t967) * t703 * t796) * t688) * t776;
t655 = (-t700 * t744 / (t789 + (t811 * t830 - t926) * pkin(3)) + (t873 * t688 - t686 + t988) * t688) * t776;
t676 = t688 * t964;
t888 = -m(3) * qJ(3,2) - mrSges(3,3);
t691 = -(t888 * pkin(2) - t1008 + t947) * t824 - (-t747 * t811 - t810 * t994 + t890) * t830;
t706 = -t730 * t830 - t878 * t824 + t845;
t753 = t843 - t888;
t786 = -qJ(3,2) / 0.2e1 + t850;
t869 = t879 * t824;
t892 = t700 * t733 * t738;
t917 = (-t708 * t807 - (t824 * t906 + (t735 * t824 + t805) * t811 + t739) * t830 + t824 * t908 - t814 ^ 2 * m(3) + (qJ(3,2) * t990) + t861) * t655 - t691 * t892 + t706 * t643 + 0.2e1 * (-t745 * t922 - (-mrSges(3,2) * t962 - t950) * t811 - t954) * t655 + t664 * t753 * t985 + (t676 * t976 + (-(t753 * pkin(2) + t1008 + t948) * t830 + ((t786 * mrSges(3,2) + t841) * t811 + (t786 * mrSges(3,1) + t842) * t810 + t763) * t976) * t703 + (t922 * t983 - 0.4e1 * t945 + (-t830 * t863 - 0.2e1 * t869) * pkin(1) + t862) * t688) * t703;
t657 = t684 + t674 / 0.2e1 + t673 / 0.2e1 - t711 / 0.2e1;
t737 = t826 * pkin(2) + pkin(3) * sin(t804);
t762 = 0.2e1 * t804;
t853 = 0.2e1 * qJ(2,1);
t800 = pkin(7) + t853;
t644 = ((cos(qJ(2,1) - pkin(7)) * t909 + cos(qJ(2,1) + t849) * t907 + t886 * t773 + t765 * t972 + t982) * t701 * t966 + (t1002 * t665 - 0.2e1 * t657 * t758 + t987 * t1004 + (-t856 * cos(t762) - t857 * cos(t853) + t871 + ((t971 - 2 * qJ(3,1)) * qJ(3,1))) * t689 / 0.2e1 + (t657 * t972 + (-cos(t800) - t811) * t689 * t910) * t1003 + (t737 + (sin(t800) * t919 + sin(t762) * t856 + sin(t853) * t857) * t966) * t704 * t797) * t689) * t777;
t656 = (-t701 * t744 / (t790 + (t811 * t832 - t925) * pkin(3)) + (t872 * t689 - t684 + t987) * t689) * t777;
t677 = t784 * t684;
t889 = -m(3) * qJ(3,1) - mrSges(3,3);
t692 = -(t889 * pkin(2) - t1009 + t947) * t826 - (-t748 * t811 - t810 * t993 + t890) * t832;
t707 = -t730 * t832 - t878 * t826 + t845;
t754 = t843 - t889;
t787 = -qJ(3,1) / 0.2e1 + t850;
t868 = t879 * t826;
t891 = t701 * t734 * t737;
t916 = (-t708 * t808 - (t826 * t906 + (t735 * t826 + t805) * t811 + t739) * t832 + t826 * t908 - t815 ^ 2 * m(3) + (qJ(3,1) * t990) + t861) * t656 - t692 * t891 + t707 * t644 + 0.2e1 * (-t745 * t921 - (-mrSges(3,2) * t961 - t950) * t811 - t954) * t656 + t665 * t754 * t984 + (t677 * t975 + (-(t754 * pkin(2) + t1009 + t948) * t832 + ((t787 * mrSges(3,2) + t841) * t811 + (t787 * mrSges(3,1) + t842) * t810 + t763) * t975) * t704 + (t921 * t983 - 0.4e1 * t944 + (-t832 * t863 - 0.2e1 * t868) * pkin(1) + t862) * t689) * t704;
t867 = t878 * t828;
t885 = -m(3) * t642 + t705 * t654 + (-(t813 * m(3) + mrSges(3,3)) * t687 / 0.2e1 + (-t730 * t822 + t867) * t702) * t986;
t866 = t878 * t830;
t884 = -m(3) * t643 + t706 * t655 + (-(t814 * m(3) + mrSges(3,3)) * t688 / 0.2e1 + (-t730 * t824 + t866) * t703) * t985;
t865 = t878 * t832;
t883 = -m(3) * t644 + t707 * t656 + (-(t815 * m(3) + mrSges(3,3)) * t689 / 0.2e1 + (-t730 * t826 + t865) * t704) * t984;
t728 = t879 * t1005 - Ifges(2,3) - Ifges(3,3) - t965;
t877 = t732 * (((t663 * t920 + t675) * t822 + (t867 + t870) * (t685 + 0.2e1 * t670 + 0.2e1 * t669 - 0.2e1 * t709) + (0.2e1 * t946 + (t708 * t822 + t970) * t828 - t1006) * t687) * t687 + t690 * t654 - t728 * t893);
t876 = t733 * (((t664 * t920 + t676) * t824 + (t866 + t869) * (t686 + 0.2e1 * t672 + 0.2e1 * t671 - 0.2e1 * t710) + (0.2e1 * t945 + (t708 * t824 + t970) * t830 - t1006) * t688) * t688 + t691 * t655 - t728 * t892);
t875 = t734 * (((t665 * t920 + t677) * t826 + (t865 + t868) * (t684 + 0.2e1 * t674 + 0.2e1 * t673 - 0.2e1 * t711) + (0.2e1 * t944 + (t708 * t826 + t970) * t832 - t1006) * t689) * t689 + t692 * t656 - t728 * t891);
t1 = [t780 * t875 + t779 * t876 + t778 * t877 - (t883 * t683 + t916 * t698) * t941 - (t884 * t681 + t917 * t697) * t942 - (t885 * t679 + t918 * t696) * t943; t783 * t875 + t782 * t876 + t781 * t877 - (t883 * t682 + t916 * t695) * t941 - (t884 * t680 + t917 * t694) * t942 - (t885 * t678 + t918 * t693) * t943; (t883 * t721 + t916 * t833) * t777 + (t884 * t720 + t917 * t831) * t776 + (t885 * t719 + t918 * t829) * t775;];
taucX  = t1;
