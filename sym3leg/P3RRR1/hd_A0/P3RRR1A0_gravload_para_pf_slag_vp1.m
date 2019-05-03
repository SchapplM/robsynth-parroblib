% Calculate Gravitation load for parallel robot
% P3RRR1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
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
%   pkin=[a2,a3,d1,d2,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:38
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRR1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(2+1,1),zeros(2+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1A0_gravload_para_pf_slag_vp1: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1A0_gravload_para_pf_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RRR1A0_gravload_para_pf_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRR1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'P3RRR1A0_gravload_para_pf_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:38:17
% EndTime: 2019-05-03 15:38:17
% DurationCPUTime: 0.33s
% Computational Cost: add. (477->90), mult. (714->168), div. (60->5), fcn. (572->20), ass. (0->74)
t671 = qJ(1,3) + qJ(2,3);
t657 = sin(t671);
t660 = cos(t671);
t677 = sin(qJ(1,3));
t680 = cos(qJ(1,3));
t703 = 0.1e1 / (t657 * t680 - t677 * t660);
t672 = qJ(1,2) + qJ(2,2);
t658 = sin(t672);
t661 = cos(t672);
t678 = sin(qJ(1,2));
t681 = cos(qJ(1,2));
t702 = 0.1e1 / (t658 * t681 - t678 * t661);
t673 = qJ(1,1) + qJ(2,1);
t659 = sin(t673);
t662 = cos(t673);
t679 = sin(qJ(1,1));
t682 = cos(qJ(1,1));
t701 = 0.1e1 / (t659 * t682 - t679 * t662);
t700 = m(2) / pkin(2);
t674 = legFrame(3,3);
t663 = sin(t674);
t666 = cos(t674);
t651 = -t663 * g(1) + t666 * g(2);
t654 = t666 * g(1) + t663 * g(2);
t615 = -t660 * (rSges(2,1) * t651 - rSges(2,2) * t654) + t657 * (rSges(2,1) * t654 + rSges(2,2) * t651);
t699 = (((-rSges(1,1) * t651 + rSges(1,2) * t654) * t680 + (rSges(1,1) * t654 + rSges(1,2) * t651) * t677) * m(1) + ((-t651 * t680 + t654 * t677) * pkin(1) + t615) * m(2)) * t703;
t675 = legFrame(2,3);
t664 = sin(t675);
t667 = cos(t675);
t652 = -t664 * g(1) + t667 * g(2);
t655 = t667 * g(1) + t664 * g(2);
t616 = -t661 * (rSges(2,1) * t652 - rSges(2,2) * t655) + t658 * (rSges(2,1) * t655 + rSges(2,2) * t652);
t698 = (((-rSges(1,1) * t652 + rSges(1,2) * t655) * t681 + (rSges(1,1) * t655 + rSges(1,2) * t652) * t678) * m(1) + ((-t652 * t681 + t655 * t678) * pkin(1) + t616) * m(2)) * t702;
t676 = legFrame(1,3);
t665 = sin(t676);
t668 = cos(t676);
t653 = -t665 * g(1) + t668 * g(2);
t656 = t668 * g(1) + t665 * g(2);
t617 = -t662 * (rSges(2,1) * t653 - rSges(2,2) * t656) + t659 * (rSges(2,1) * t656 + rSges(2,2) * t653);
t697 = (((-rSges(1,1) * t653 + rSges(1,2) * t656) * t682 + (rSges(1,1) * t656 + rSges(1,2) * t653) * t679) * m(1) + ((-t653 * t682 + t656 * t679) * pkin(1) + t617) * m(2)) * t701;
t696 = t615 * t703;
t695 = t616 * t702;
t694 = t617 * t701;
t630 = t666 * t657 + t663 * t660;
t631 = -t657 * t663 + t666 * t660;
t632 = t667 * t658 + t664 * t661;
t633 = -t658 * t664 + t667 * t661;
t634 = t668 * t659 + t665 * t662;
t635 = -t659 * t665 + t668 * t662;
t693 = 0.1e1 / pkin(1);
t691 = koppelP(1,1);
t690 = koppelP(2,1);
t689 = koppelP(3,1);
t688 = koppelP(1,2);
t687 = koppelP(2,2);
t686 = koppelP(3,2);
t685 = rSges(3,1);
t684 = rSges(3,2);
t683 = xP(3);
t670 = cos(t683);
t669 = sin(t683);
t650 = -t669 * t688 + t670 * t691;
t649 = -t669 * t687 + t670 * t690;
t648 = -t669 * t686 + t670 * t689;
t647 = -t669 * t691 - t670 * t688;
t646 = -t669 * t690 - t670 * t687;
t645 = -t669 * t689 - t670 * t686;
t623 = pkin(1) * (-t679 * t665 + t668 * t682) + t635 * pkin(2);
t622 = pkin(1) * (-t678 * t664 + t667 * t681) + t633 * pkin(2);
t621 = pkin(1) * (-t677 * t663 + t666 * t680) + t631 * pkin(2);
t620 = pkin(1) * (t665 * t682 + t668 * t679) + t634 * pkin(2);
t619 = pkin(1) * (t664 * t681 + t667 * t678) + t632 * pkin(2);
t618 = pkin(1) * (t663 * t680 + t666 * t677) + t630 * pkin(2);
t1 = [-m(3) * g(1) + (t631 * t699 + t633 * t698 + t635 * t697 + (-t621 * t696 - t622 * t695 - t623 * t694) * t700) * t693; -m(3) * g(2) + (t630 * t699 + t632 * t698 + t634 * t697 + (-t618 * t696 - t619 * t695 - t620 * t694) * t700) * t693; ((g(1) * t685 + g(2) * t684) * t669 + (g(1) * t684 - g(2) * t685) * t670) * m(3) + ((t634 * t650 + t635 * t647) * t697 + (t632 * t649 + t633 * t646) * t698 + (t630 * t648 + t631 * t645) * t699 + (-(t620 * t650 + t623 * t647) * t694 - (t619 * t649 + t622 * t646) * t695 - (t618 * t648 + t621 * t645) * t696) * t700) * t693;];
taugX  = t1;
