Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace '.stl trim','' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'upper','upper'.ToUpper() }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'lower','lower'.ToUpper() }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'maxillary','upper'.ToUpper() }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'mandibular','lower'.ToUpper() }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'retainer','RETAINER' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'refinement','REFINEMENT' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'ö','o' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'ä','a' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'ü','u' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'ß','s' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace '_','' }
Get-ChildItem -File | Rename-Item -NewName { 
   ForEach($Name in $_.Name) 
	{
     		If ($Name -match "LOWER")  
     			{$Name -replace " Step ", "LOWER"}
      		elseif($Name -match "UPPER") 
      			{$Name -replace " Step ", "UPPER" }
		elseif($Name -match "MAXILLARY") 
      			{$Name -replace " Step ", "UPPER" }
		elseif($Name -match "MANDIBULAR") 
      			{$Name -replace " Step ", "LOWER" }
		elseif($Name -match "RETAINER") 
      			{$Name -replace " Step ", "_RETAINER_" }
		elseif($Name -match "REFINEMENT") 
      			{$Name -replace " Step ", "_REFINEMENT_" } 
	}
   
}
Get-ChildItem -exclude *.ps1 | Rename-Item -NewName {
ForEach($string1 in $_.Name) 
	{
			$splitstring = $string1 -split 'of'
			$first = $splitstring[0]
			$last = $splitstring[1]
     		If ($last -match "pts")  
     			{$string1 -replace $last, ".pts"}
		elseif($last -match "stl") {$string1 -replace $last, ".stl"}
                  
      		 
	}
}

Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'of','' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace ' ','' }


Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER0.PTS','LN00.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER0.PTS','UN00.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER1.PTS','LN01.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER1.PTS','UN01.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER2.PTS','LN02.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER2.PTS','UN02.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER3.PTS','LN03.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER3.PTS','UN03.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER4.PTS','LN04.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER4.PTS','UN04.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER5.PTS','LN05.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER5.PTS','UN05.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER6.PTS','LN06.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER6.PTS','UN06.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER7.PTS','LN07.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER7.PTS','UN07.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER8.PTS','LN08.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER8.PTS','UN08.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER9.PTS','LN09.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER9.PTS','UN09.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER00.PTS','LN00.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER00.PTS','UN00.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER01.PTS','LN01.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER01.PTS','UN01.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER02.PTS','LN02.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER02.PTS','UN02.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER03.PTS','LN03.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER03.PTS','UN03.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER04.PTS','LN04.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER04.PTS','UN04.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER05.PTS','LN05.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER05.PTS','UN05.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER06.PTS','LN06.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER06.PTS','UN06.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER07.PTS','LN07.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER07.PTS','UN07.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER08.PTS','LN08.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER08.PTS','UN08.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER09.PTS','LN09.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER09.PTS','UN09.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER10.PTS','LN10.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER10.PTS','UN10.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER11.PTS','LN11.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER11.PTS','UN11.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER12.PTS','LN12.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER12.PTS','UN12.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER13.PTS','LN13.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER13.PTS','UN13.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER14.PTS','LN14.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER14.PTS','UN14.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER15.PTS','LN15.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER15.PTS','UN15.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER16.PTS','LN16.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER16.PTS','UN16.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER17.PTS','LN17.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER17.PTS','UN17.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER18.PTS','LN18.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER18.PTS','UN18.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER19.PTS','LN19.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER19.PTS','UN19.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER20.PTS','LN20.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER20.PTS','UN20.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER21.PTS','LN21.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER21.PTS','UN21.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER22.PTS','LN22.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER22.PTS','UN22.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER23.PTS','LN23.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER23.PTS','UN23.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER24.PTS','LN24.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER24.PTS','UN24.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER25.PTS','LN25.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER25.PTS','UN25.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER26.PTS','LN26.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER26.PTS','UN26.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER27.PTS','LN27.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER27.PTS','UN27.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER28.PTS','LN28.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER28.PTS','UN28.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER29.PTS','LN29.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER29.PTS','UN29.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER30.PTS','LN30.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER30.PTS','UN30.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER31.PTS','LN31.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER31.PTS','UN31.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER32.PTS','LN32.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER32.PTS','UN32.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER33.PTS','LN33.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER33.PTS','UN33.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER34.PTS','LN34.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER34.PTS','UN34.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER35.PTS','LN35.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER35.PTS','UN35.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER36.PTS','LN36.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER36.PTS','UN36.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER37.PTS','LN37.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER37.PTS','UN37.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER38.PTS','LN38.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER38.PTS','UN38.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER39.PTS','LN39.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER39.PTS','UN39.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER40.PTS','LN40.PTS' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER40.PTS','UN40.PTS' }


Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER0.stl','LN00.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER0.stl','UN00.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER1.stl','LN01.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER1.stl','UN01.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER2.stl','LN02.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER2.stl','UN02.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER3.stl','LN03.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER3.stl','UN03.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER4.stl','LN04.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER4.stl','UN04.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER5.stl','LN05.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER5.stl','UN05.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER6.stl','LN06.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER6.stl','UN06.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER7.stl','LN07.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER7.stl','UN07.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER8.stl','LN08.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER8.stl','UN08.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER9.stl','LN09.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER9.stl','UN09.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER00.stl','LN00.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER00.stl','UN00.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER01.stl','LN01.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER10.stl','UN01.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER02.stl','LN02.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER02.stl','UN02.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER03.stl','LN03.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER03.stl','UN03.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER04.stl','LN04.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER04.stl','UN04.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER05.stl','LN05.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER05.stl','UN05.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER06.stl','LN06.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER06.stl','UN06.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER07.stl','LN07.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER07.stl','UN07.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER08.stl','LN08.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER08.stl','UN08.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER09.stl','LN09.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER09.stl','UN09.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER10.stl','LN10.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER10.stl','UN10.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER11.stl','LN11.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER11.stl','UN11.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER12.stl','LN12.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER12.stl','UN12.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER13.stl','LN13.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER13.stl','UN13.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER14.stl','LN14.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER14.stl','UN14.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER15.stl','LN15.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER15.stl','UN15.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER16.stl','LN16.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER16.stl','UN16.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER17.stl','LN17.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER17.stl','UN17.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER18.stl','LN18.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER18.stl','UN18.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER19.stl','LN19.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER19.stl','UN19.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER20.stl','LN20.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER20.stl','UN20.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER21.stl','LN21.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER21.stl','UN21.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER22.stl','LN22.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER22.stl','UN22.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER23.stl','LN23.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER23.stl','UN23.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER24.stl','LN24.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER24.stl','UN24.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER25.stl','LN25.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER25.stl','UN25.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER26.stl','LN26.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER26.stl','UN26.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER27.stl','LN27.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER27.stl','UN27.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER28.stl','LN28.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER28.stl','UN28.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER29.stl','LN29.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER29.stl','UN29.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER30.stl','LN30.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER30.stl','UN30.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER31.stl','LN31.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER31.stl','UN31.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER32.stl','LN32.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER32.stl','UN32.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER33.stl','LN33.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER33.stl','UN33.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER34.stl','LN34.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER34.stl','UN34.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER35.stl','LN35.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER35.stl','UN35.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER36.stl','LN36.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER36.stl','UN36.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER37.stl','LN37.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER37.stl','UN37.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER38.stl','LN38.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER38.stl','UN38.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER39.stl','LN39.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER39.stl','UN39.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER40.stl','LN40.stl' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER40.stl','UN40.stl' }

Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER0','LN00' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER0','UN00' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER1','LN01' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER1','UN01' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER2','LN02' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER2','UN02' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER3','LN03' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER3','UN03' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER4','LN04' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER4','UN04' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER5','LN05' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER5','UN05' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER6','LN06' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER6','UN06' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER7','LN07' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER7','UN07' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER8','LN08' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER8','UN08' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER9','LN09' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER9','UN09' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER00','LN00' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER00','UN00' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER01','LN01' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER01','UN01' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER02','LN02' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER02','UN02' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER03','LN03' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER03','UN03' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER04','LN04' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER04','UN04' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER05','LN05' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER05','UN05' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER06','LN06' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER06','UN06' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER07','LN07' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER07','UN07' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER08','LN08' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER08','UN08' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER09','LN09' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER09','UN09' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER10','LN10' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER10','UN10' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER11','LN11' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER11','UN11' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER12','LN12' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER12','UN12' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER13','LN13' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER13','UN13' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER14','LN14' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER14','UN14' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER15','LN15' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER15','UN15' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER16','LN16' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER16','UN16' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER17','LN17' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER17','UN17' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER18','LN18' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER18','UN18' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER19','LN19' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER19','UN19' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER20','LN20' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER20','UN20' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER21','LN21' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER21','UN21' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER22','LN22' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER22','UN22' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER23','LN23' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER23','UN23' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER24','LN24' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER24','UN24' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER25','LN25' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER25','UN25' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER26','LN26' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER26','UN26' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER27','LN27' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER27','UN27' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER28','LN28' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER28','UN28' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER29','LN29' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER29','UN29' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER30','LN30' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER30','UN30' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER31','LN31' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER31','UN31' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER32','LN32' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER32','UN32' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER33','LN33' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER33','UN33' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER34','LN34' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER34','UN34' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER35','LN35' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER35','UN35' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER36','LN36' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER36','UN36' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER37','LN37' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER37','UN37' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER38','LN38' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER38','UN38' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER39','LN39' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER39','UN39' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'LOWER40','LN40' }
Get-ChildItem -File | Rename-Item -NewName { $_.Name -replace 'UPPER40','UN40' }
