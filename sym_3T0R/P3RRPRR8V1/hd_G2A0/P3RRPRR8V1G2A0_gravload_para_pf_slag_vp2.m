% Calculate Gravitation load for parallel robot
% P3RRPRR8V1G2A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-04 17:05
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:04:29
% EndTime: 2022-11-04 17:04:29
% DurationCPUTime: 0.43s
% Computational Cost: add. (555->115), mult. (789->190), div. (30->9), fcn. (600->23), ass. (0->86)
t817 = sin(pkin(5));
t874 = pkin(2) * t817;
t822 = legFrame(3,2);
t808 = sin(t822);
t811 = cos(t822);
t794 = t811 * g(1) - t808 * g(2);
t819 = pkin(4) + qJ(3,3);
t814 = 0.1e1 / t819;
t826 = sin(qJ(1,3));
t832 = cos(qJ(1,3));
t873 = (t826 * g(3) - t832 * t794) * t814;
t823 = legFrame(2,2);
t809 = sin(t823);
t812 = cos(t823);
t795 = t812 * g(1) - t809 * g(2);
t820 = pkin(4) + qJ(3,2);
t815 = 0.1e1 / t820;
t828 = sin(qJ(1,2));
t834 = cos(qJ(1,2));
t872 = (t828 * g(3) - t834 * t795) * t815;
t824 = legFrame(1,2);
t810 = sin(t824);
t813 = cos(t824);
t796 = t813 * g(1) - t810 * g(2);
t821 = pkin(4) + qJ(3,1);
t816 = 0.1e1 / t821;
t830 = sin(qJ(1,1));
t836 = cos(qJ(1,1));
t871 = (t830 * g(3) - t836 * t796) * t816;
t818 = cos(pkin(5));
t790 = m(3) * pkin(1) + mrSges(3,1) * t818 - mrSges(3,2) * t817 + mrSges(2,1);
t791 = t808 * g(1) + t811 * g(2);
t800 = t817 * mrSges(3,1) + t818 * mrSges(3,2) + mrSges(2,2);
t825 = sin(qJ(2,3));
t831 = cos(qJ(2,3));
t845 = g(3) * t832 + t794 * t826;
t848 = t831 * pkin(1) + pkin(2) * cos(qJ(2,3) + pkin(5));
t870 = 0.1e1 / t848 * ((t845 * t790 + t791 * t800) * t825 - (t791 * t790 - t845 * t800) * t831);
t792 = t809 * g(1) + t812 * g(2);
t827 = sin(qJ(2,2));
t833 = cos(qJ(2,2));
t844 = g(3) * t834 + t795 * t828;
t847 = t833 * pkin(1) + pkin(2) * cos(qJ(2,2) + pkin(5));
t869 = 0.1e1 / t847 * ((t844 * t790 + t792 * t800) * t827 - (t792 * t790 - t844 * t800) * t833);
t793 = t810 * g(1) + t813 * g(2);
t829 = sin(qJ(2,1));
t835 = cos(qJ(2,1));
t843 = g(3) * t836 + t796 * t830;
t846 = t835 * pkin(1) + pkin(2) * cos(qJ(2,1) + pkin(5));
t868 = 0.1e1 / t846 * ((t843 * t790 + t793 * t800) * t829 - (t793 * t790 - t843 * t800) * t835);
t804 = t818 * pkin(2) + pkin(1);
t867 = t804 * t811;
t866 = t804 * t812;
t865 = t804 * t813;
t864 = t808 * t804;
t863 = t809 * t804;
t862 = t810 * t804;
t858 = mrSges(2,3) + mrSges(3,3) - mrSges(1,2);
t801 = m(3) * qJ(3,3) + t858;
t839 = -t790 * t831 + t800 * t825 - mrSges(1,1);
t861 = t814 * ((-t801 * t826 + t839 * t832) * t794 + (-t801 * t832 - t839 * t826) * g(3));
t802 = m(3) * qJ(3,2) + t858;
t838 = -t790 * t833 + t800 * t827 - mrSges(1,1);
t860 = t815 * ((-t802 * t828 + t838 * t834) * t795 + (-t802 * t834 - t838 * t828) * g(3));
t803 = m(3) * qJ(3,1) + t858;
t837 = -t790 * t835 + t800 * t829 - mrSges(1,1);
t859 = t816 * ((-t803 * t830 + t837 * t836) * t796 + (-t803 * t836 - t837 * t830) * g(3));
t857 = t811 * t874;
t856 = t812 * t874;
t855 = t813 * t874;
t854 = t808 * t874;
t853 = t809 * t874;
t852 = t810 * t874;
t842 = t804 * t831 - t825 * t874;
t851 = 0.1e1 / t842 * t861;
t841 = t804 * t833 - t827 * t874;
t850 = 0.1e1 / t841 * t860;
t840 = t804 * t835 - t829 * t874;
t849 = 0.1e1 / t840 * t859;
t789 = t829 * t804 + t835 * t874;
t788 = t827 * t804 + t833 * t874;
t787 = t825 * t804 + t831 * t874;
t780 = -t836 * t821 + t840 * t830;
t779 = -t834 * t820 + t841 * t828;
t778 = -t832 * t819 + t842 * t826;
t1 = [((t830 * t865 + t852) * t835 + (-t830 * t855 + t862) * t829) * t849 + t810 * t868 + ((t828 * t866 + t853) * t833 + (-t828 * t856 + t863) * t827) * t850 + t809 * t869 + ((t826 * t867 + t854) * t831 + (-t826 * t857 + t864) * t825) * t851 + t808 * t870 - g(1) * m(4) + (-(t780 * t813 + t810 * t789) * t871 - (t779 * t812 + t809 * t788) * t872 - (t778 * t811 + t808 * t787) * t873) * m(3); ((-t830 * t862 + t855) * t835 + (t830 * t852 + t865) * t829) * t849 + t813 * t868 + ((-t828 * t863 + t856) * t833 + (t828 * t853 + t866) * t827) * t850 + t812 * t869 + ((-t826 * t864 + t857) * t831 + (t826 * t854 + t867) * t825) * t851 + t811 * t870 - g(2) * m(4) + (-(-t780 * t810 + t813 * t789) * t871 - (-t779 * t809 + t812 * t788) * t872 - (-t778 * t808 + t811 * t787) * t873) * m(3); t832 * t861 + t834 * t860 + t836 * t859 - g(3) * m(4) + (-(t830 * t821 + t846 * t836) * t871 - (t828 * t820 + t847 * t834) * t872 - (t826 * t819 + t848 * t832) * t873) * m(3);];
taugX  = t1;
