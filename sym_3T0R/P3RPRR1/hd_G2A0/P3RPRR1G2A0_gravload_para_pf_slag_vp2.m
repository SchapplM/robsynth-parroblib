% Calculate Gravitation load for parallel robot
% P3RPRR1G2A0
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
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRR1G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:24:57
% EndTime: 2020-03-09 21:24:57
% DurationCPUTime: 0.30s
% Computational Cost: add. (579->132), mult. (564->133), div. (27->4), fcn. (330->66), ass. (0->87)
t647 = m(2) + m(3);
t636 = legFrame(3,2);
t626 = sin(t636);
t629 = cos(t636);
t578 = t629 * g(1) - t626 * g(2);
t655 = pkin(7) + qJ(3,3);
t617 = qJ(1,3) + t655;
t610 = cos(t617);
t643 = mrSges(3,2) * g(3);
t645 = mrSges(3,1) * g(3);
t566 = -t610 * (mrSges(3,1) * t578 - t643) + sin(t617) * (mrSges(3,2) * t578 + t645);
t581 = 0.1e1 / (pkin(1) * sin(t655) + sin(qJ(3,3)) * pkin(2));
t673 = t566 * t581;
t637 = legFrame(2,2);
t627 = sin(t637);
t630 = cos(t637);
t579 = t630 * g(1) - t627 * g(2);
t656 = pkin(7) + qJ(3,2);
t618 = qJ(1,2) + t656;
t611 = cos(t618);
t567 = -t611 * (mrSges(3,1) * t579 - t643) + sin(t618) * (mrSges(3,2) * t579 + t645);
t582 = 0.1e1 / (pkin(1) * sin(t656) + sin(qJ(3,2)) * pkin(2));
t672 = t567 * t582;
t638 = legFrame(1,2);
t628 = sin(t638);
t631 = cos(t638);
t580 = t631 * g(1) - t628 * g(2);
t657 = pkin(7) + qJ(3,1);
t619 = qJ(1,1) + t657;
t612 = cos(t619);
t568 = -t612 * (mrSges(3,1) * t580 - t643) + sin(t619) * (mrSges(3,2) * t580 + t645);
t583 = 0.1e1 / (pkin(1) * sin(t657) + sin(qJ(3,1)) * pkin(2));
t671 = t568 * t583;
t670 = (t626 * g(1) + t629 * g(2)) * t647;
t669 = (t627 * g(1) + t630 * g(2)) * t647;
t668 = (t628 * g(1) + t631 * g(2)) * t647;
t603 = t647 * pkin(1) + mrSges(1,1);
t596 = t603 * g(3);
t632 = m(3) * pkin(2) + mrSges(2,1);
t613 = t632 * g(3);
t633 = qJ(1,3) + pkin(7);
t614 = cos(t633);
t639 = cos(qJ(1,3));
t642 = g(3) * mrSges(2,2);
t644 = mrSges(1,2) * g(3);
t667 = t581 * ((-t578 * t632 + t642) * t614 + (t578 * mrSges(2,2) + t613) * sin(t633) + (-t578 * t603 + t644) * t639 + sin(qJ(1,3)) * (t578 * mrSges(1,2) + t596) + t566);
t634 = qJ(1,2) + pkin(7);
t615 = cos(t634);
t640 = cos(qJ(1,2));
t666 = t582 * ((-t579 * t632 + t642) * t615 + (t579 * mrSges(2,2) + t613) * sin(t634) + (-t579 * t603 + t644) * t640 + sin(qJ(1,2)) * (t579 * mrSges(1,2) + t596) + t567);
t635 = qJ(1,1) + pkin(7);
t616 = cos(t635);
t641 = cos(qJ(1,1));
t665 = t583 * ((-t580 * t632 + t642) * t616 + (t580 * mrSges(2,2) + t613) * sin(t635) + (-t580 * t603 + t644) * t641 + sin(qJ(1,1)) * (t580 * mrSges(1,2) + t596) + t568);
t648 = 0.1e1 / pkin(3);
t664 = t648 / 0.2e1;
t604 = t636 + t633;
t597 = qJ(3,3) + t604;
t605 = -t636 + t633;
t598 = qJ(3,3) + t605;
t663 = sin(t597) + sin(t598);
t606 = t637 + t634;
t599 = qJ(3,2) + t606;
t607 = -t637 + t634;
t600 = qJ(3,2) + t607;
t662 = sin(t599) + sin(t600);
t608 = t638 + t635;
t601 = qJ(3,1) + t608;
t609 = -t638 + t635;
t602 = qJ(3,1) + t609;
t661 = sin(t601) + sin(t602);
t660 = -cos(t598) + cos(t597);
t659 = -cos(t600) + cos(t599);
t658 = -cos(t602) + cos(t601);
t654 = t667 / 0.2e1;
t653 = t666 / 0.2e1;
t652 = t665 / 0.2e1;
t651 = t664 * t673;
t650 = t664 * t672;
t649 = t664 * t671;
t625 = qJ(1,1) - t638;
t624 = qJ(1,1) + t638;
t623 = qJ(1,2) - t637;
t622 = qJ(1,2) + t637;
t621 = qJ(1,3) - t636;
t620 = qJ(1,3) + t636;
t1 = [t661 * t652 - t628 * t668 - (t661 * pkin(3) + (sin(t608) + sin(t609)) * pkin(2) + (sin(t624) + sin(t625)) * pkin(1)) * t649 + t662 * t653 - t627 * t669 - (t662 * pkin(3) + (sin(t606) + sin(t607)) * pkin(2) + (sin(t622) + sin(t623)) * pkin(1)) * t650 + t663 * t654 - t626 * t670 - (t663 * pkin(3) + (sin(t604) + sin(t605)) * pkin(2) + (sin(t620) + sin(t621)) * pkin(1)) * t651 - g(1) * m(4); t658 * t652 - t631 * t668 + (-t658 * pkin(3) + (cos(t609) - cos(t608)) * pkin(2) + (cos(t625) - cos(t624)) * pkin(1)) * t649 + t659 * t653 - t630 * t669 + (-t659 * pkin(3) + (cos(t607) - cos(t606)) * pkin(2) + (cos(t623) - cos(t622)) * pkin(1)) * t650 + t660 * t654 - t629 * t670 + (-t660 * pkin(3) + (cos(t605) - cos(t604)) * pkin(2) + (cos(t621) - cos(t620)) * pkin(1)) * t651 - g(2) * m(4); t610 * t667 + t611 * t666 + t612 * t665 - g(3) * m(4) + ((-t641 * pkin(1) - pkin(2) * t616 - pkin(3) * t612) * t671 + (-t640 * pkin(1) - pkin(2) * t615 - pkin(3) * t611) * t672 + (-t639 * pkin(1) - pkin(2) * t614 - pkin(3) * t610) * t673) * t648;];
taugX  = t1;
