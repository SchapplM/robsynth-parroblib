% Calculate Gravitation load for parallel robot
% P3RPRR1G3A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRR1G3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRR1G3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:26:48
% EndTime: 2020-03-09 21:26:48
% DurationCPUTime: 0.33s
% Computational Cost: add. (579->132), mult. (663->140), div. (27->4), fcn. (330->66), ass. (0->89)
t728 = m(2) + m(3);
t761 = m(1) * rSges(1,2);
t760 = m(2) * rSges(2,2);
t759 = rSges(3,1) * g(3);
t758 = rSges(3,2) * g(3);
t757 = m(3) / pkin(3);
t721 = legFrame(3,2);
t712 = sin(t721);
t715 = cos(t721);
t660 = t715 * g(1) - t712 * g(2);
t737 = pkin(7) + qJ(3,3);
t703 = qJ(1,3) + t737;
t693 = sin(t703);
t648 = t693 * (rSges(3,1) * t660 - t758) + cos(t703) * (rSges(3,2) * t660 + t759);
t663 = 0.1e1 / (pkin(1) * sin(t737) + sin(qJ(3,3)) * pkin(2));
t756 = t648 * t663;
t722 = legFrame(2,2);
t713 = sin(t722);
t716 = cos(t722);
t661 = t716 * g(1) - t713 * g(2);
t738 = pkin(7) + qJ(3,2);
t704 = qJ(1,2) + t738;
t694 = sin(t704);
t649 = t694 * (rSges(3,1) * t661 - t758) + cos(t704) * (rSges(3,2) * t661 + t759);
t664 = 0.1e1 / (pkin(1) * sin(t738) + sin(qJ(3,2)) * pkin(2));
t755 = t649 * t664;
t723 = legFrame(1,2);
t714 = sin(t723);
t717 = cos(t723);
t662 = t717 * g(1) - t714 * g(2);
t739 = pkin(7) + qJ(3,1);
t705 = qJ(1,1) + t739;
t695 = sin(t705);
t650 = t695 * (rSges(3,1) * t662 - t758) + cos(t705) * (rSges(3,2) * t662 + t759);
t665 = 0.1e1 / (pkin(1) * sin(t739) + sin(qJ(3,1)) * pkin(2));
t754 = t650 * t665;
t753 = (t712 * g(1) + t715 * g(2)) * t728;
t752 = (t713 * g(1) + t716 * g(2)) * t728;
t751 = (t714 * g(1) + t717 * g(2)) * t728;
t679 = m(1) * rSges(1,1) + t728 * pkin(1);
t666 = g(3) * t679;
t699 = m(2) * rSges(2,1) + m(3) * pkin(2);
t686 = g(3) * t699;
t718 = qJ(1,3) + pkin(7);
t700 = sin(t718);
t724 = sin(qJ(1,3));
t740 = g(3) * t760;
t741 = g(3) * t761;
t750 = t663 * ((t660 * t760 + t686) * cos(t718) + (t699 * t660 - t740) * t700 + (t660 * t761 + t666) * cos(qJ(1,3)) + (t679 * t660 - t741) * t724 + t648 * m(3));
t719 = qJ(1,2) + pkin(7);
t701 = sin(t719);
t725 = sin(qJ(1,2));
t749 = t664 * ((t661 * t760 + t686) * cos(t719) + (t699 * t661 - t740) * t701 + (t661 * t761 + t666) * cos(qJ(1,2)) + (t679 * t661 - t741) * t725 + t649 * m(3));
t720 = qJ(1,1) + pkin(7);
t702 = sin(t720);
t726 = sin(qJ(1,1));
t748 = t665 * ((t662 * t760 + t686) * cos(t720) + (t699 * t662 - t740) * t702 + (t662 * t761 + t666) * cos(qJ(1,1)) + (t679 * t662 - t741) * t726 + t650 * m(3));
t687 = t721 + t718;
t680 = qJ(3,3) + t687;
t688 = -t721 + t718;
t681 = qJ(3,3) + t688;
t747 = -sin(t680) + sin(t681);
t689 = t722 + t719;
t682 = qJ(3,2) + t689;
t690 = -t722 + t719;
t683 = qJ(3,2) + t690;
t746 = -sin(t682) + sin(t683);
t691 = t723 + t720;
t684 = qJ(3,1) + t691;
t692 = -t723 + t720;
t685 = qJ(3,1) + t692;
t745 = -sin(t684) + sin(t685);
t744 = cos(t681) + cos(t680);
t743 = cos(t683) + cos(t682);
t742 = cos(t685) + cos(t684);
t736 = t757 / 0.2e1;
t735 = t750 / 0.2e1;
t734 = t749 / 0.2e1;
t733 = t748 / 0.2e1;
t732 = t736 * t756;
t731 = t736 * t755;
t730 = t736 * t754;
t711 = qJ(1,1) - t723;
t710 = qJ(1,1) + t723;
t709 = qJ(1,2) - t722;
t708 = qJ(1,2) + t722;
t707 = qJ(1,3) - t721;
t706 = qJ(1,3) + t721;
t1 = [t742 * t733 - t714 * t751 + (-t742 * pkin(3) + (-cos(t692) - cos(t691)) * pkin(2) + (-cos(t710) - cos(t711)) * pkin(1)) * t730 + t743 * t734 - t713 * t752 + (-t743 * pkin(3) + (-cos(t690) - cos(t689)) * pkin(2) + (-cos(t709) - cos(t708)) * pkin(1)) * t731 + t744 * t735 - t712 * t753 + (-t744 * pkin(3) + (-cos(t688) - cos(t687)) * pkin(2) + (-cos(t707) - cos(t706)) * pkin(1)) * t732 - m(4) * g(1); t745 * t733 - t717 * t751 - (t745 * pkin(3) + (-sin(t691) + sin(t692)) * pkin(2) + (-sin(t710) + sin(t711)) * pkin(1)) * t730 + t746 * t734 - t716 * t752 - (t746 * pkin(3) + (-sin(t689) + sin(t690)) * pkin(2) + (-sin(t708) + sin(t709)) * pkin(1)) * t731 + t747 * t735 - t715 * t753 - (t747 * pkin(3) + (-sin(t687) + sin(t688)) * pkin(2) + (-sin(t706) + sin(t707)) * pkin(1)) * t732 - m(4) * g(2); -t693 * t750 - t694 * t749 - t695 * t748 - m(4) * g(3) + ((t726 * pkin(1) + pkin(2) * t702 + pkin(3) * t695) * t754 + (t725 * pkin(1) + pkin(2) * t701 + pkin(3) * t694) * t755 + (t724 * pkin(1) + pkin(2) * t700 + pkin(3) * t693) * t756) * t757;];
taugX  = t1;
