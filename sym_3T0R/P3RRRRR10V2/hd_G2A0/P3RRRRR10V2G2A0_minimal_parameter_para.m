% Return the minimum parameter vector for
% P3RRRRR10V2G2A0
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% 
% Output:
% MPV [18x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 02:06
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = P3RRRRR10V2G2A0_minimal_parameter_para(pkin, m, mrSges, Ifges, koppelP)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6),zeros(3,3)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3RRRRR10V2G2A0_minimal_parameter_para: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR10V2G2A0_minimal_parameter_para: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR10V2G2A0_minimal_parameter_para: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR10V2G2A0_minimal_parameter_para: Ifges has to be [4x6] (double)'); 
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V2G2A0_minimal_parameter_para: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From minimal_parameter_parrob_matlab.m
t5720 = m(2) + m(3);
t5725 = (t5720 * pkin(6));
t5724 = 2 * pkin(7) * mrSges(3,3) + Ifges(3,2);
t5723 = m(3) * pkin(7) + mrSges(3,3);
t5722 = (pkin(2) ^ 2);
t5721 = pkin(7) ^ 2;
t5719 = sin(pkin(4));
t1 = [pkin(1) ^ 2 * t5720 + Ifges(1,3) + (m(3) * t5722 + Ifges(2,2) + (2 * mrSges(2,3) + t5725) * pkin(6)) * t5719 ^ 2; pkin(1) * t5720 + mrSges(1,1); mrSges(1,2) + (-mrSges(2,3) - t5725) * t5719; Ifges(2,1) - Ifges(2,2) + (t5721 - t5722) * m(3) + t5724; t5723 * pkin(2) + Ifges(2,4); Ifges(2,5); Ifges(2,6); Ifges(2,3) + (t5721 + t5722) * m(3) + t5724; m(3) * pkin(2) + mrSges(2,1); mrSges(2,2) - t5723; Ifges(3,1) - Ifges(3,2); Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3); mrSges(3,1); mrSges(3,2); m(4);];
MPV  = t1;
