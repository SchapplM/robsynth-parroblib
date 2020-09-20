% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRR1G2P2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x8]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRR1G2P2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:30
% EndTime: 2020-03-09 21:18:30
% DurationCPUTime: 0.30s
% Computational Cost: add. (865->82), mult. (468->156), div. (126->5), fcn. (558->18), ass. (0->78)
t789 = legFrame(3,2);
t780 = sin(t789);
t786 = pkin(7) + qJ(2,3);
t777 = qJ(3,3) + t786;
t762 = sin(t777);
t765 = cos(t777);
t771 = sin(t786);
t774 = cos(t786);
t814 = 0.1e1 / (-t774 * t762 + t765 * t771);
t820 = t780 * t814;
t790 = legFrame(2,2);
t781 = sin(t790);
t787 = pkin(7) + qJ(2,2);
t778 = qJ(3,2) + t787;
t763 = sin(t778);
t766 = cos(t778);
t772 = sin(t787);
t775 = cos(t787);
t813 = 0.1e1 / (-t775 * t763 + t766 * t772);
t819 = t781 * t813;
t791 = legFrame(1,2);
t782 = sin(t791);
t788 = pkin(7) + qJ(2,1);
t779 = qJ(3,1) + t788;
t764 = sin(t779);
t767 = cos(t779);
t773 = sin(t788);
t776 = cos(t788);
t812 = 0.1e1 / (-t776 * t764 + t767 * t773);
t818 = t782 * t812;
t783 = cos(t789);
t817 = t783 * t814;
t784 = cos(t790);
t816 = t784 * t813;
t785 = cos(t791);
t815 = t785 * t812;
t811 = t814 * (-pkin(2) * t774 - pkin(3) * t765);
t810 = t814 * t765;
t809 = t813 * (-pkin(2) * t775 - pkin(3) * t766);
t808 = t813 * t766;
t807 = t812 * (-pkin(2) * t776 - pkin(3) * t767);
t806 = t812 * t767;
t750 = pkin(2) * t771 + pkin(3) * t762;
t805 = t750 * t817;
t804 = t762 * t820;
t751 = pkin(2) * t772 + pkin(3) * t763;
t803 = t751 * t816;
t802 = t763 * t819;
t752 = pkin(2) * t773 + pkin(3) * t764;
t801 = t752 * t815;
t800 = t764 * t818;
t799 = t750 * t820;
t798 = t762 * t817;
t797 = t751 * t819;
t796 = t763 * t816;
t795 = t752 * t818;
t794 = t764 * t815;
t793 = 0.1e1 / pkin(2);
t792 = 0.1e1 / pkin(3);
t761 = -t782 * g(1) - t785 * g(2);
t760 = -t785 * g(1) + t782 * g(2);
t759 = -t781 * g(1) - t784 * g(2);
t758 = -t784 * g(1) + t781 * g(2);
t757 = -t780 * g(1) - t783 * g(2);
t756 = -t783 * g(1) + t780 * g(2);
t740 = g(3) * t776 - t760 * t773;
t739 = g(3) * t775 - t758 * t772;
t738 = g(3) * t774 - t756 * t771;
t737 = g(3) * t773 + t760 * t776;
t736 = g(3) * t772 + t758 * t775;
t735 = g(3) * t771 + t756 * t774;
t734 = g(3) * t767 - t760 * t764;
t733 = g(3) * t766 - t758 * t763;
t732 = g(3) * t765 - t756 * t762;
t731 = g(3) * t764 + t760 * t767;
t730 = g(3) * t763 + t758 * t766;
t729 = g(3) * t762 + t756 * t765;
t1 = [t780 * t757 + t781 * t759 + t782 * t761, 0, (-t735 * t798 - t736 * t796 - t737 * t794) * t793, (-t738 * t798 - t739 * t796 - t740 * t794) * t793, 0, (-t729 * t798 - t730 * t796 - t731 * t794 + (t729 * t805 + t730 * t803 + t731 * t801) * t792) * t793, (-t732 * t798 - t733 * t796 - t734 * t794 + (t732 * t805 + t733 * t803 + t734 * t801) * t792) * t793, -g(1); t783 * t757 + t784 * t759 + t785 * t761, 0, (t735 * t804 + t736 * t802 + t737 * t800) * t793, (t738 * t804 + t739 * t802 + t740 * t800) * t793, 0, (t729 * t804 + t730 * t802 + t731 * t800 + (-t729 * t799 - t730 * t797 - t731 * t795) * t792) * t793, (t732 * t804 + t733 * t802 + t734 * t800 + (-t732 * t799 - t733 * t797 - t734 * t795) * t792) * t793, -g(2); 0, 0, (-t735 * t810 - t736 * t808 - t737 * t806) * t793, (-t738 * t810 - t739 * t808 - t740 * t806) * t793, 0, (-t729 * t810 - t730 * t808 - t731 * t806 + (-t729 * t811 - t730 * t809 - t731 * t807) * t792) * t793, (-t732 * t810 - t733 * t808 - t734 * t806 + (-t732 * t811 - t733 * t809 - t734 * t807) * t792) * t793, -g(3);];
tau_reg  = t1;
